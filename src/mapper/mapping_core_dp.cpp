// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.0
// Date    : 2017-12-24
// License : GNU GPL 3
// *******************************************************************************************


#include "mapping_core.h"
#include "../libs/asmlib.h"
#include "../common/utils.h"

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) ((i) & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

//#define LEV_MYERS_ASSERTIONS

//**********************************************************************************************************
//
//**********************************************************************************************************
uint32_t CMappingCore::search_in_malicious_group_with_DP(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
	uint32_t* sub_SA, uint32_t sub_sa_size,
	CAux_arrays* aux_arr, CAux_params* aux_params, bool* aux_struct_valid,
	bool dir_gen_flag, uint32_t max_mismatches,
	uint32_t sub_SA_index_offs)
{
	if (sensitive_mode)
		return max_mismatches + 1;

	vector<uint32_t> *sum_of_localizations;
	vector<uchar_t> *list_merger_tmp;

	LevMyers* levMyers;
	levMyers = (pattern_len < 256)
		? (pattern_len < 128) ? levMyers128 : levMyers256
		: levMyers64;

	if (dir_gen_flag)
	{
		sum_of_localizations = &sum_of_localizations_dir;
		list_merger_tmp = &list_merger_tmp_dir;
	}
	else
	{
		sum_of_localizations = &sum_of_localizations_rc;
		list_merger_tmp = &list_merger_tmp_rc;
	}

	if (!(*aux_struct_valid))
	{
		build_aux_structures_for_DP(sub_SA, sub_sa_size, min_read_len, aux_arr, aux_params, dir_gen_flag);
		*aux_struct_valid = true;
		sum_of_localizations->resize(sub_sa_size);
		list_merger_tmp->assign(sub_sa_size + 1, 0);
	}

	//-----------------------------------------------------------------------------------------------
	//searching starts
	//-----------------------------------------------------------------------------------------------
	uint32_t cur_match = 0;
	uint32_t best_match = max_mismatches + 1;

	if (dir_gen_flag)
		for (uint32_t k = 0; k < no_of_substrings; k++)
		{
			cur_substr_values_from_pattern[k] = get_substring_value_from_pattern(aux_params->substring_length, pattern_ptr,
				(*(aux_params->substr_pos_arr))[k]);
		}
	else
		for (uint32_t k = 0; k < no_of_substrings; k++)
		{
			//the third parameter is the starting point of a substring for reversed pattern
			cur_substr_values_from_pattern[k] = get_substring_value_from_pattern(aux_params->substring_length, pattern_ptr,
				pattern_len - (*(aux_params->substr_pos_arr))[k] - aux_params->substring_length);
		}

	in_list_desc.clear();

	for (uint32_t k_DP = 0; k_DP < aux_params->no_of_substrings_for_DP; k_DP++)
	{
		uint32_t aux_index_read = (*(aux_params->substr_pos_mapping_arr))[k_DP].first;
		uint32_t aux_index_ref = (*(aux_params->substr_pos_mapping_arr))[k_DP].second;
		if (aux_index_ref >= aux_params->no_of_substrings_in_text_for_DP)
		{
			cerr << "\n" << aux_index_read << "  " << aux_index_ref << "   " << aux_params->no_of_substrings_in_text_for_DP << "  " << k_DP <<
				"  " << aux_params->no_of_substrings_for_DP << "\n";
		}

		uint32_t where_is_index = aux_arr[aux_index_ref].whereToPut_array[cur_substr_values_from_pattern[aux_index_read]];

		in_list_desc.push_back(make_pair(aux_arr[aux_index_ref].localization_array.begin() + where_is_index,
			aux_arr[aux_index_ref].localization_array.begin() + where_is_index + aux_arr[aux_index_ref].howManySubstr_array[cur_substr_values_from_pattern[aux_index_read]]));
	}

	// Filtering out too long lists
	uint32_t req_occ = no_of_substrings - stage_major;
	sort(in_list_desc.begin(), in_list_desc.end(),
		[](const pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator> &x, const pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator> &y)
	{return distance(x.first, x.second) < distance(y.first, y.second);});

	if (req_occ > 3)
	{
		uint32_t sum = 0;
		for (auto p = in_list_desc.begin(); p != in_list_desc.end() - 1; ++p)
			sum += (uint32_t)(p->second - p->first);
		if (in_list_desc.back().second - in_list_desc.back().first > sum)
		{
			in_list_desc.pop_back();
			--req_occ;
		}
	}

	vector<uint32_t>::iterator q = list_merger(in_list_desc, req_occ, list_merger_tmp, (*sum_of_localizations).begin());

	if (max_mismatches >= lev_alg_thr)
	{
#ifdef LEV_MYERS_ASSERTIONS
		edit_dist->LevMyers_PP(pattern_ptr, pattern_len, genome_t::direct);		// Direction is always direct as read is already reverse-complemented (if dir_falg is false)
#endif
		levMyers->preprocess(pattern_ptr, pattern_len, genome_t::direct);
	}
	for (auto p = (*sum_of_localizations).begin(); p != q; ++p)
	{
		//-------------------------------DIR GENOME ----------------------------------------
		if (dir_gen_flag)
		{
			if (identical_sequences_dir[*p] == false)
			{
				if (test_CF_128b_SSE_difference(CF_value_SSE_for_a_dir_pattern, CF_vector_128b_SSE_dir_gen[*p + sub_SA_index_offs],
					4 * max_mismatches))
				{
					ref_pos_t ref_pos;
					if (max_mismatches >= lev_alg_thr) {
						levMyers->dynamicProgramming(sub_SA[*p] - start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos, cur_match);
#ifdef LEV_MYERS_ASSERTIONS
						uint32_t ref_pos2, cur_match2;
						edit_dist->LevMyers(sub_SA[*p] - start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos2, cur_match2);
						ASSERT((ref_pos == ref_pos2) && (cur_match == cur_match2), "CMappingCore::search_in_malicious_group_with_DP() error!");
#endif
					}
					else
						cur_match = edit_dist->LevDiag(pattern_ptr, pattern_len, sub_SA[*p], start_pos, end_pos, max_mismatches, 1);

					no_of_CF_accepted_in_malicious_group++;
					if (cur_match <= max_mismatches)
						no_of_Lev_positive++;
				}
				else
				{
					cur_match = max_mismatches + 1;
					no_of_CF_discarded_in_malicious_group++;
				}
			}
		}
		//--------------------------------RC GENOME ---------------------------------------
		else
		{
			if (identical_sequences_rc[*p] == false)
			{
				if (test_CF_128b_SSE_difference(CF_value_SSE_for_a_rc_pattern, CF_vector_128b_SSE_rc_gen[*p + sub_SA_index_offs], 4 * max_mismatches))
				{
					ref_pos_t ref_pos;
					if (max_mismatches >= lev_alg_thr) {
						levMyers->dynamicProgramming(ref_size - sub_SA[*p] - pattern_len + start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos, cur_match);
#ifdef LEV_MYERS_ASSERTIONS
						uint32_t ref_pos2, cur_match2;
						edit_dist->LevMyers(ref_size - sub_SA[*p] - pattern_len + start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos2, cur_match2);
						ASSERT(ref_pos == ref_pos2 && cur_match == cur_match2, "CMappingCore::search_in_malicious_group_with_DP() error!");
#endif
					}
					else
						cur_match = edit_dist->LevDiag(pattern_ptr, pattern_len, sub_SA[*p], pattern_len - end_pos, pattern_len - start_pos, max_mismatches, 0);

					no_of_CF_accepted_in_malicious_group++;
					if (cur_match <= max_mismatches)
						no_of_Lev_positive++;
				}
				else
				{
					cur_match = max_mismatches + 1;
					no_of_CF_discarded_in_malicious_group++;
				}
			}
		}
		//--------------------------------------------------------------------------------

		if (cur_match < best_match)
			best_match = cur_match;

		if (cur_match <= max_mismatches)
		{
			if (dir_gen_flag)
				dir_match_positions.push_back(make_pair(sub_SA[*p] - CUR_OFFSET, cur_match));
			else
				rc_match_positions.push_back(make_pair(ref_size - (sub_SA[*p] - CUR_OFFSET) - 1, cur_match));
		}
	}

	return best_match;
}

//**********************************************************************************************************
//
//**********************************************************************************************************
void CMappingCore::build_aux_structures_for_DP(uint32_t* sub_SA, uint32_t sub_size, uint32_t pattern_len, CAux_arrays* aux_arr, CAux_params* aux_params, bool dir_gen_flag)
{
	uint32_t substring_value;
	uint32_t start_substring_position;

	for (uint32_t i = 1; i < group_size_substr_length_mapping.size(); ++i)
		if (group_size_substr_length_mapping[i] > sub_size)
		{
			aux_params->substring_length = i - 1;
			aux_params->substr_pos_arr = &(substr_pos_all_lengths[aux_params->substring_length]);
			aux_params->substr_pos_in_text_arr = &(substr_pos_in_text_all_lengths[aux_params->substring_length]);
			aux_params->substr_pos_mapping_arr = &(substr_pos_mapping_all_lengths[aux_params->substring_length]);
			aux_params->no_of_substrings_for_DP = (uint32_t)aux_params->substr_pos_mapping_arr->size();
			aux_params->no_of_substrings_in_text_for_DP = (uint32_t)aux_params->substr_pos_in_text_arr->size();
			break;
		}

	aux_params->LEN = 1 << 2 * aux_params->substring_length;

	for (uint32_t i = 0; i < aux_params->no_of_substrings_in_text_for_DP; i++)
	{
		aux_arr[i].howManySubstr_array.assign(aux_params->LEN + 1, 0);
		aux_arr[i].whereToPut_array.assign(aux_params->LEN + 1, 0);
		aux_arr[i].localization_array.assign(sub_size, 0);
	}

	uint32_t* substr_values_aux_arr = new uint32_t[sub_size * aux_params->no_of_substrings_in_text_for_DP];
	A_memset(substr_values_aux_arr, 0, sub_size * aux_params->no_of_substrings_in_text_for_DP * sizeof(uint32_t));

	old_pos_in_text_dir = 0;
	old_pos_in_text_rc = 0;

	if (dir_gen_flag)
		identical_sequences_dir.assign(sub_size, 0);
	else
		identical_sequences_rc.assign(sub_size, 0);


	for (uint32_t i = 0; i < sub_size; i++) //count the occurrences of individual substring_values
	{
		//------------------------------------------------------------------------------------------------
		//
		//------------------------------------------------------------------------------------------------

		uint32_t left_end_pos;
#if 0
		bool cmp;
#endif

		if (dir_gen_flag)
		{
			left_end_pos = sub_SA[i] - CUR_OFFSET;
#if 0
			cmp = compare_str_to_str_in_text(left_end_pos - stage_major, pattern_len + 2 * stage_major, old_pos_in_text_dir);
#endif
		}
		else
		{
			left_end_pos = ref_size - sub_SA[i] - pattern_len + CUR_OFFSET;
#if 0
			cmp = compare_str_to_str_in_text(left_end_pos - stage_major, pattern_len + 2 * stage_major, old_pos_in_text_rc);
#endif
		}

		uint32_t val_offset = i*aux_params->no_of_substrings_in_text_for_DP;

#if 0
		// Future: possible optimization. If the strings are equal, we can use the substrings from the previous one instead of calculating them
		if(!cmp)
		{
#endif
			//-------------------------------------------------------------------------------------------
			//sequences are different. We have to copy left_end_pos i.e. the starting point of the new seq.
			//-------------------------------------------------------------------------------------------
			if (dir_gen_flag)
			{
				old_pos_in_text_dir = left_end_pos - stage_major;
				uint32_t sa_pos = sub_SA[i] - CUR_OFFSET;
				for (uint32_t k_DP = 0; k_DP < aux_params->no_of_substrings_in_text_for_DP; k_DP++)
				{
					start_substring_position = (*(aux_params->substr_pos_in_text_arr))[k_DP];
					substring_value = get_substring_value(aux_params->substring_length, (int32_t)start_substring_position, sa_pos);

					++aux_arr[k_DP].howManySubstr_array[substring_value];
					substr_values_aux_arr[val_offset + k_DP] = substring_value;
				}
			}
			else
			{
				old_pos_in_text_rc = left_end_pos - stage_major;
				uint32_t sa_pos = ref_size - sub_SA[i] - pattern_len + CUR_OFFSET;
				for (uint32_t k_DP = 0; k_DP < aux_params->no_of_substrings_in_text_for_DP; k_DP++)
				{
					start_substring_position = pattern_len - (*(aux_params->substr_pos_in_text_arr))[k_DP] - aux_params->substring_length;
					substring_value = get_substring_value(aux_params->substring_length, (int32_t)start_substring_position, sa_pos);

					++aux_arr[k_DP].howManySubstr_array[substring_value];
					substr_values_aux_arr[val_offset + k_DP] = substring_value;
				}
			}
#if 0
		}
		else	//the sequence is the same as a previous one. Substring values can be copied, no need to evaluate them
		{
			if (dir_gen_flag)
				identical_sequences_dir[i] = true;
			else
				identical_sequences_rc[i] = true;

			for (uint32_t k_DP = 0; k_DP < aux_params->no_of_substrings_in_text_for_DP; k_DP++)
			{
				//				substr_values_aux_arr[i*no_of_substrings_for_DP + k_DP] = substr_values_aux_arr[(i-1)*no_of_substrings_for_DP + k_DP];
				//				substr_values_aux_arr[i*aux_params->no_of_substrings_in_text_for_DP + k_DP] = substr_values_aux_arr[(i-1)*aux_params->no_of_substrings_in_text_for_DP + k_DP];
				substr_values_aux_arr[val_offset + k_DP] = substr_values_aux_arr[val_offset - aux_params->no_of_substrings_in_text_for_DP + k_DP];
				//				++aux_arr[k_DP].howManySubstr_array[substr_values_aux_arr[(i-1)*no_of_substrings_for_DP + k_DP]];
				//				++aux_arr[k_DP].howManySubstr_array[substr_values_aux_arr[(i-1)*aux_params->no_of_substrings_in_text_for_DP + k_DP]];
				++aux_arr[k_DP].howManySubstr_array[substr_values_aux_arr[val_offset - aux_params->no_of_substrings_in_text_for_DP + k_DP]];
				//--------------++aux_arr[k].howManySubstr_array[substring_value];
				//--------------substr_values_aux_arr[i][k] = substring_value;
			}
		}
#endif
	} // end_for(uint32_t i = 0; i < sub_size; i++) 

	for (uint32_t k_DP = 0; k_DP < aux_params->no_of_substrings_in_text_for_DP; k_DP++)
	{
		aux_arr[k_DP].whereToPut_array[0] = aux_arr[k_DP].howManySubstr_array[0];
		for (uint32_t i = 1; i <= aux_params->LEN; i++)
			aux_arr[k_DP].whereToPut_array[i] = aux_arr[k_DP].howManySubstr_array[i] + aux_arr[k_DP].whereToPut_array[i - 1];

		uint32_t val_offset = (sub_size - 1) * aux_params->no_of_substrings_in_text_for_DP + k_DP;
		for (uint32_t i = sub_size; i > 0; i--, val_offset -= aux_params->no_of_substrings_in_text_for_DP)
		{
			aux_arr[k_DP].localization_array[aux_arr[k_DP].whereToPut_array[substr_values_aux_arr[val_offset]] - 1] = i - 1;
			aux_arr[k_DP].whereToPut_array[substr_values_aux_arr[val_offset]]--;
		}

		aux_arr[k_DP].howManySubstr_array[aux_params->LEN] = 0;		//Search for illegal substring_values does not take place
	}

	delete[] substr_values_aux_arr;
}

//**********************************************************************************************************
void CMappingCore::calculate_CF_128b_SEE_for_dir_and_rc_gen(uint32_t fixedLeft, uint32_t fixedRight, uint32_t pattern_len, uint32_t position_read_from_SA_array, Vec16c& CF, bool dir_gen_flag)
{
	uint32_t left_end;

	if (dir_gen_flag)
		left_end = position_read_from_SA_array - CUR_OFFSET;
	else
		left_end = ref_size - position_read_from_SA_array - pattern_len + CUR_OFFSET;

	CF = 0;
	uchar_t* T_ptr = ref_ptr + (left_end >> 1);

	if (fixedLeft)
	{
		uint32_t len_in_chars = fixedLeft;

		if ((left_end & 1) == 0)
		{
			T_ptr = ref_ptr + (left_end >> 1);
		}
		else
		{
			T_ptr = ref_ptr + (left_end >> 1);
			uchar_t byte = LO_NIBBLE(*T_ptr) << 4;
			T_ptr++;
			byte += HI_NIBBLE(*T_ptr);

			CF = add_saturated(CF, CF_SSE_lut[byte]);
			--len_in_chars;
		}

		uint32_t no_full_bytes = (len_in_chars - 1) >> 1;

		for (uint32_t i = 0; i < no_full_bytes; ++i)
		{
			CF = add_saturated(CF, CF_SSE_lut[*T_ptr]);
			uchar_t byte = LO_NIBBLE(*T_ptr) << 4;
			T_ptr++;
			byte += HI_NIBBLE(*T_ptr);
			CF = add_saturated(CF, CF_SSE_lut[byte]);
		}

		if ((len_in_chars & 1) == 0)
			CF = add_saturated(CF, CF_SSE_lut[*T_ptr]);
	}
	//-------------------------------
	//outside matching part:
	//-------------------------------
	if (fixedRight != pattern_len)
	{
		uint32_t len_in_chars = pattern_len - fixedRight;

		left_end += fixedRight;

		if ((left_end & 1) == 0)
		{
			T_ptr = ref_ptr + (left_end >> 1);
		}
		else
		{
			T_ptr = ref_ptr + (left_end >> 1);
			uchar_t byte = LO_NIBBLE(*T_ptr) << 4;
			T_ptr++;
			byte += HI_NIBBLE(*T_ptr);

			CF = add_saturated(CF, CF_SSE_lut[byte]);
			--len_in_chars;
		}

		uint32_t no_full_bytes = (len_in_chars - 1) >> 1;

		for (uint32_t i = 0; i < no_full_bytes; ++i)
		{
			CF = add_saturated(CF, CF_SSE_lut[*T_ptr]);
			uchar_t byte = LO_NIBBLE(*T_ptr) << 4;
			T_ptr++;
			byte += HI_NIBBLE(*T_ptr);
			CF = add_saturated(CF, CF_SSE_lut[byte]);
		}

		if ((len_in_chars & 1) == 0)
			CF = add_saturated(CF, CF_SSE_lut[*T_ptr]);
	}
}

//**********************************************************************************************************
void CMappingCore::calculate_CF_128b_SSE_for_a_pattern(uchar_t *pattern, uchar_t* pattern_sft, uint32_t pattern_len, uint32_t fixedLeft, uint32_t fixedRight, Vec16c& CF_value_SSE)
{
	uchar_t *p, *q;

	CF_value_SSE = 0;

	//---------------------CF before matching part:
	if (fixedLeft)
	{
		uint32_t len_in_chars = fixedLeft;
		uint32_t no_full_bytes = (len_in_chars - 1) / 2;

		p = &pattern[0];
		q = &pattern_sft[1];

		for (uint32_t i = 0; i < no_full_bytes; ++i)
		{
			CF_value_SSE = add_saturated(CF_value_SSE, CF_SSE_lut[*p++]);
			CF_value_SSE = add_saturated(CF_value_SSE, CF_SSE_lut[*q++]);
		}

		if ((len_in_chars & 1) == 0)
			CF_value_SSE = add_saturated(CF_value_SSE, CF_SSE_lut[*p++]);
	}

	//--------------------outside matching part
	if (fixedRight != pattern_len)
	{
		uint32_t len_in_chars = pattern_len - fixedRight;
		uint32_t no_full_bytes = (len_in_chars - 1) / 2;
		uint32_t fixedRight_div2 = fixedRight / 2;

		if (fixedRight & 1)
		{
			p = &pattern_sft[fixedRight_div2 + 1];
			q = &pattern[fixedRight_div2 + 1];
		}
		else
		{
			p = &pattern[fixedRight_div2];
			q = &pattern_sft[fixedRight_div2 + 1];
		}

		for (uint32_t i = 0; i < no_full_bytes; ++i)
		{
			CF_value_SSE = add_saturated(CF_value_SSE, CF_SSE_lut[*p++]);
			CF_value_SSE = add_saturated(CF_value_SSE, CF_SSE_lut[*q++]);
		}

		if ((len_in_chars & 1) == 0)
			CF_value_SSE = add_saturated(CF_value_SSE, CF_SSE_lut[*p++]);
	}
}

//**********************************************************************************************************
bool CMappingCore::test_CF_128b_SSE_difference(Vec16c CF_1, Vec16c CF_2, uint32_t max_diff)
{
	uint32_t CF_diff = 0;

	CF_diff = horizontal_add(abs(CF_1 - CF_2));

	return CF_diff <= max_diff;
}

// EOF 