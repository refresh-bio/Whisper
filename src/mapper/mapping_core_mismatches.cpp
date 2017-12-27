#include "mapping_core.h"
#include <algorithm>
#include <utility>
#include <math.h>
#include "../libs/asmlib.h"
#include "../asm_common/utils.h"

#ifdef MISMATCHES_CODE
//*****************************************************************************************
uint32_t CMappingCore::search_in_malicious_group_dir_and_rc(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
									   uint32_t* sub_SA, uint32_t sub_sa_size, 
									   CAux_arrays* aux_arr, CAux_params* aux_params, bool* aux_struct_valid, bool dir_gen_flag, uint32_t max_mismatches,
									   uint32_t sub_SA_index_offs)
{
	vector<uint32_t> *sum_of_localizations;
	vector<uchar_t> *list_merger_tmp;
	
	if(dir_gen_flag)
	{
		sum_of_localizations = &sum_of_localizations_dir;
		list_merger_tmp = &list_merger_tmp_dir;
	}
	else
	{
		sum_of_localizations = &sum_of_localizations_rc;
		list_merger_tmp = &list_merger_tmp_rc;
	}

	if(!(*aux_struct_valid))
	{

		//build_aux_structures_dir_and_rc(sub_SA, sub_sa_size, pattern_len, aux_arr, aux_params, dir_gen_flag);
		//POPRAWIAM: 29.09.2014
		build_aux_structures_dir_and_rc(sub_SA, sub_sa_size, min_read_len, aux_arr, aux_params, dir_gen_flag);
		*aux_struct_valid = true;
		sum_of_localizations->resize(sub_sa_size);
//		list_merger_tmp->resize(sub_sa_size+1);
//		A_memset(list_merger_tmp->data(), 0, list_merger_tmp->size());
		list_merger_tmp->assign(sub_sa_size+1, 0);
	}
	//-----------------------------------------------------------------------------------------------
	//Teraz bedzie szukanie
	//-----------------------------------------------------------------------------------------------
	uint32_t cur_match = 0;
	uint32_t best_match = max_mismatches + 1;
//	int32_t cmp = 0;
//	uint32_t count = 0;

	if(dir_gen_flag)
		for(uint32_t k = 0; k < no_of_substrings; k++)
		{
			cur_substr_values_from_pattern[k] = get_substring_value_from_pattern(aux_params->substring_length, pattern_ptr,   
				(*(aux_params->substr_pos_arr))[k]);
//			count += aux_arr[k].howManySubstr_array[cur_substr_values_from_pattern[k]];
		}
	else
		for(uint32_t k = 0; k < no_of_substrings; k++)
		{
			//the third parameter is the starting point of a substring for reversed pattern
			cur_substr_values_from_pattern[k] = get_substring_value_from_pattern(aux_params->substring_length, pattern_ptr,  
				pattern_len - (*(aux_params->substr_pos_arr))[k] - aux_params->substring_length);
//			count += aux_arr[k].howManySubstr_array[cur_substr_values_from_pattern[k]];
		}

//	if(count > 0.00 * sub_sa_size)		// tymczasowe w³¹czenie, ze zawsze sie robi merge
	{
		in_list_desc.clear();
		for(uint32_t k = 0; k < no_of_substrings; k++)
		{
			uint32_t where_is_index = aux_arr[k].whereToPut_array[cur_substr_values_from_pattern[k]];
			in_list_desc.push_back(make_pair(aux_arr[k].localization_array.begin() + where_is_index,
				aux_arr[k].localization_array.begin() + where_is_index + aux_arr[k].howManySubstr_array[cur_substr_values_from_pattern[k]]));
		}
//-------
		// Filtering out too long lists
		sort(in_list_desc.begin(), in_list_desc.end(), 
			[](const pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator> &x, const pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator> &y)
				{return distance(x.first, x.second) < distance(y.first, y.second);});
//				{return x.second-x.first < y.second-y.first;});
		while(in_list_desc.size() > stage_major+3)
		{
			uint32_t sum = 0;
			for(auto p = in_list_desc.begin(); p != in_list_desc.end()-1; ++p)
				sum += (uint32_t) (p->second - p->first);
			if(in_list_desc.back().second - in_list_desc.back().first > sum)
				in_list_desc.pop_back();
			else
				break;
		}
//-----------
		vector<uint32_t>::iterator q = list_merger(in_list_desc, (uint32_t) (in_list_desc.size()-stage_major), list_merger_tmp, (*sum_of_localizations).begin());
		for(auto p = (*sum_of_localizations).begin(); p != q; ++p)
		{
			//-------------------------------DIR GENOME ----------------------------------------
			if(dir_gen_flag)
			{
				if(identical_sequences_dir[*p] == false)

						cur_match = test_selected_data(pattern_ptr, pattern_sft_ptr, pattern_len, sub_SA, *p, max_mismatches);

			}

			//--------------------------------RC GENOME ---------------------------------------
			else
			{
				if(identical_sequences_rc[*p] == false)

					cur_match = test_selected_rc_data(pattern_ptr, pattern_sft_ptr, pattern_len, sub_SA, *p, max_mismatches);

			}
			//--------------------------------------------------------------------------------

			if(cur_match < best_match)
				best_match = cur_match;

			if(cur_match <= max_mismatches)
			{
				if(dir_gen_flag)
					dir_match_positions.push_back(make_pair(sub_SA[*p] - CUR_OFFSET, cur_match));
				else
					rc_match_positions.push_back(make_pair(ref_size - (sub_SA[*p] - CUR_OFFSET) - 1, cur_match));
			}
		}
	}
//	else
//	{
////		for(uint32_t k = 0; k < no_of_substrings; k++)
//		// Tu szukac po listach, ktore maja najmniej elementów
//		for(uint32_t k = 0; k <= stage_major; k++)
//		{
//			uint32_t where_is_index = aux_arr[k].whereToPut_array[cur_substr_values_from_pattern[k]];//to powinno siê odbywaæ warunkowo, je¿eli wartoœæ wystepuje
//			
//			for(uint32_t i = 0; i < aux_arr[k].howManySubstr_array[cur_substr_values_from_pattern[k]]; i++)
//			{
//				uint32_t where_is_in_SA = aux_arr[k].localization_array[where_is_index];
//				where_is_index++;	
//
//				//-------------------------------DIR GENOME ----------------------------------------
//				if(dir_gen_flag)
//				{
//					if(identical_sequences_dir[where_is_in_SA] == false)
//#ifdef USE_32b_CF_AND_MISMATCHES
//					{
//						if(test_CF_difference(CF_value_for_a_dir_pattern, CF_array_dir_gen[where_is_in_SA + sub_SA_index_offs], 2 * max_mismatches))
//						{
//							cur_match = test_selected_data(pattern_ptr, pattern_sft_ptr, pattern_len, sub_SA, where_is_in_SA, stage_major);
//							no_of_CF_accepted_in_malicious_group++;
//						}
//						else
//						{
//							cur_match = max_mismatches + 1;
//							no_of_CF_discarded_in_malicious_group++;
//						}
//					}
//#else
//							cur_match = test_selected_data(pattern_ptr, pattern_sft_ptr, pattern_len, sub_SA, where_is_in_SA, stage_major);
//#endif
//				}
//
//				//--------------------------------RC GENOME ---------------------------------------
//				else
//				{
//					if(identical_sequences_rc[where_is_in_SA] == false)
//#ifdef USE_32b_CF_AND_MISMATCHES
//					{
//						if(test_CF_difference(CF_value_for_a_rc_pattern, CF_array_rc_gen[where_is_in_SA + sub_SA_index_offs], 2 * max_mismatches))
//							cur_match = test_selected_rc_data(pattern_ptr, pattern_sft_ptr, pattern_len, sub_SA, where_is_in_SA, stage_major);
//						else
//							cur_match = max_mismatches;
//					}
//#else
//							cur_match = test_selected_rc_data(pattern_ptr, pattern_sft_ptr, pattern_len, sub_SA, where_is_in_SA, stage_major);
//#endif
//				}
	//			//--------------------------------------------------------------------------------

	//			if(cur_match < best_match)
	//				best_match = cur_match;

	//			if(cur_match <= stage_major)
	//			{
	//				if(dir_gen_flag)
	//					dir_match_positions.push_back(make_pair(sub_SA[where_is_in_SA] - CUR_OFFSET, cur_match));
	//				else
	//					rc_match_positions.push_back(make_pair(ref_size - (sub_SA[where_is_in_SA] + CUR_OFFSET) - 1, cur_match));
	//			}
	//		}
	//	}
	//}
	
	return best_match;
};

//***********************************************************************

//**************************************************************************************************
void CMappingCore::build_aux_structures_dir_and_rc(uint32_t* sub_SA, uint32_t sub_size, uint32_t pattern_len, CAux_arrays* aux_arr, CAux_params* aux_params, bool dir_gen_flag)
{
	uint32_t substring_value;
	uint32_t start_substring_position; 

	//!!!!UWAGA! PRZEROBIÆ AUX_ARR!!!
//----------------------------------------

	for(uint32_t i = 1; i < group_size_substr_length_mapping.size(); ++i)
		if(group_size_substr_length_mapping[i] > sub_size)
		{
			aux_params->substring_length = i-1;
			aux_params->substr_pos_arr = &(substr_pos_all_lengths[aux_params->substring_length]);
			break;
		}
		
	aux_params->LEN = 1 << 2*aux_params->substring_length;

//-----------------------------------------

	for(uint32_t i = 0; i < no_of_substrings; i++)
	{
		aux_arr[i].howManySubstr_array.assign (aux_params->LEN + 1, 0);
		aux_arr[i].whereToPut_array.assign (aux_params->LEN + 1, 0);
		aux_arr[i].localization_array.assign (sub_size, 0);		
	}
	
	uint32_t* substr_values_aux_arr = new uint32_t[sub_size * no_of_substrings];		
	A_memset(substr_values_aux_arr, 0, sub_size * no_of_substrings * sizeof(uint32_t));

	old_pos_in_text_dir = 0;
	old_pos_in_text_rc = 0;


	if(dir_gen_flag)
		identical_sequences_dir.assign(sub_size, 0);
	else
		identical_sequences_rc.assign(sub_size, 0);


	for(uint32_t i = 0; i < sub_size; i++) //zliczanie wyst¹pieñ poszczególnych substring_values
	{
		//------------------------------------------------------------------------------------------------
		//
		//------------------------------------------------------------------------------------------------

		uint32_t left_end_pos;
		bool cmp;
		uint32_t val_offset = i*no_of_substrings;

		if(dir_gen_flag)
		{
			left_end_pos = sub_SA[i] - CUR_OFFSET;
			//zmieniæ na max_read_len:
			cmp = compare_str_to_str_in_text(left_end_pos, pattern_len, old_pos_in_text_dir);
		}
		else
		{
			left_end_pos = ref_size - sub_SA[i] - pattern_len + CUR_OFFSET;
			//zmieniæ na max_read_len:
			cmp = compare_str_to_str_in_text(left_end_pos, pattern_len, old_pos_in_text_rc);
		}

		//!!!!!!!!!!!!!!!!!!!!!!!!POPRAWIAM NA CHWILÊ, SZUKAJAC BLEDU, 30.09.2014
		//if(!cmp)
		if(1)
		{
			//-------------------------------------------------------------------------------------------
			//sequences are different. We have to copy left_end_pos i.e. the starting point if the new seq.
			//-------------------------------------------------------------------------------------------
			if(dir_gen_flag)
			{
				old_pos_in_text_dir = left_end_pos;
				uint32_t sa_pos = sub_SA[i] - CUR_OFFSET;
				for(uint32_t k = 0; k < no_of_substrings; k++)
				{
					start_substring_position = (*(aux_params->substr_pos_arr))[k];
					substring_value = get_substring_value(aux_params->substring_length, start_substring_position, sa_pos);

					++aux_arr[k].howManySubstr_array[substring_value];
					substr_values_aux_arr[val_offset + k] = substring_value;		
				}
			} //end_for(uint32 k = 0; k < no_of_substrings; k++)
			else
			{
				old_pos_in_text_rc = left_end_pos;
				uint32_t sa_pos = ref_size - sub_SA[i] - pattern_len + CUR_OFFSET;
				for(uint32_t k = 0; k < no_of_substrings; k++)
				{
				   start_substring_position = pattern_len - (*(aux_params->substr_pos_arr))[k] - aux_params->substring_length;
				   substring_value = get_substring_value(aux_params->substring_length, start_substring_position, sa_pos);
	
					++aux_arr[k].howManySubstr_array[substring_value];
					substr_values_aux_arr[val_offset + k] = substring_value;		
				}
			} //end_for(uint32 k = 0; k < no_of_substrings; k++)

/*			if(dir_gen_flag)
				old_pos_in_text_dir = left_end_pos;
			else
				old_pos_in_text_rc = left_end_pos;

			for(uint32_t k = 0; k < no_of_substrings; k++)
			{
				if(dir_gen_flag)
				{
				   start_substring_position = (*(aux_params->substr_pos_arr))[k];
				   substring_value = get_substring_value(aux_params->substring_length, start_substring_position, sub_SA[i] - CUR_OFFSET);
				}
				else
				{
				   start_substring_position = pattern_len - (*(aux_params->substr_pos_arr))[k] - aux_params->substring_length;
				   substring_value = get_substring_value(aux_params->substring_length, start_substring_position, ref_size - sub_SA[i] - pattern_len + CUR_OFFSET);
				}
	
				++aux_arr[k].howManySubstr_array[substring_value];
				substr_values_aux_arr[val_offset + k] = substring_value;		
			} //end_for(uint32 k = 0; k < no_of_substrings; k++)
*/
		}
		else	//the sequence is the same as a previous one. Substring values can be copied, no need to evaluate them
		{
			if(dir_gen_flag)
				identical_sequences_dir[i] = true;
			else
				identical_sequences_rc[i] = true;

			for(uint32_t k = 0; k < no_of_substrings; k++)
			{
//				substr_values_aux_arr[i*no_of_substrings + k] = substr_values_aux_arr[(i-1)*no_of_substrings + k];
				substr_values_aux_arr[val_offset + k] = substr_values_aux_arr[val_offset - no_of_substrings + k];
//				++aux_arr[k].howManySubstr_array[substr_values_aux_arr[(i-1)*no_of_substrings + k]];
				++aux_arr[k].howManySubstr_array[substr_values_aux_arr[val_offset - no_of_substrings + k]];
				//--------------++aux_arr[k].howManySubstr_array[substring_value];
				//--------------substr_values_aux_arr[i][k] = substring_value;
			}
		}
	} // end_for(uint32_t i = 0; i < sub_size; i++) 

	//------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------
	for(uint32_t k = 0; k < no_of_substrings; k++ )
	{
		aux_arr[k].whereToPut_array[0] = aux_arr[k].howManySubstr_array[0];
		for(uint32_t i = 1; i <= aux_params->LEN; i++)
			aux_arr[k].whereToPut_array[i] = aux_arr[k].howManySubstr_array[i] + aux_arr[k].whereToPut_array[i-1];
	
		for(uint32_t i = sub_size; i > 0; i--)
		{
			aux_arr[k].localization_array[aux_arr[k].whereToPut_array[substr_values_aux_arr[(i - 1)*no_of_substrings + k]] - 1] = i - 1;
			aux_arr[k].whereToPut_array[substr_values_aux_arr[(i - 1)*no_of_substrings + k]]--;
		}

		aux_arr[k].howManySubstr_array[aux_params->LEN] = 0;		//nie chcemy, ¿eby odbywa³o sie szukanie dla niedozwolonych wartoœci substring_values
	}
		
	delete [] substr_values_aux_arr;
}

//*******************************************************************************************
//**************************************************************************************
uint32_t CMappingCore::test_selected_data(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len, uint32_t* sa_part, uint32_t SA_index, uint32_t max_mismatches)
{//-------------------------------ONLY DIR GENOME --------------------------------------------------

//	uint32_t length_before_matching_part = CUR_OFFSET;
	//ma tak zostac dla RZD:
	uint32_t length_after_matching_part  = pattern_len - predef_match_length - CUR_OFFSET;
	uint32_t cur_match = 0;
	uint32_t done_mismatches = 0;


	//-------------------------------------------------------------------------------------------------
	//check left part of the read, on the left from a matching part
	//-----------------------------------------------------------------------------------------------
	if(CUR_OFFSET)
	{
		uint32_t left_index = 0;					//left index of the current segment
		uint32_t left_offset = CUR_OFFSET - left_index;

		//one mismatch per segment is allowed for all segments placed on the right of the current one.
		//That gives: (stage_minor - 1)* 1 mismatches
		uint32_t allowed_mismatches_for_the_current_segment = max_mismatches - (stage_minor - 1); 
																								
		for(uint32_t i = 0; i < stage_minor; i++)
		{
			uint32_t right_index = segment_segments[i].second;		//right index of the current segment
			uint32_t length = right_index - left_index;
			uint32_t left_end_index_in_dir_gen = sa_part[SA_index] - left_offset;

			switch(left_end_index_in_dir_gen & 1)
			{
			case 0:
				if(left_index & 1)	//parzysta pozycja w tekscie, nieparzysty offset
					cur_match = check_how_many_mismatches(	ref_ptr + (left_end_index_in_dir_gen>>1),
																	pattern_sft_ptr + ((left_index + 1)>>1),					
																	length, false, allowed_mismatches_for_the_current_segment);
				else                //parzysta pozycja w tekscie, parzysty offset
					cur_match = check_how_many_mismatches(	ref_ptr + (left_end_index_in_dir_gen>>1),
																	pattern_ptr + (left_index >> 1),					
																	length, false, allowed_mismatches_for_the_current_segment);
				break;
			case 1:
				if(left_index & 1)	//nieparzysta pozycja w tekscie, nieparzysty offset
					cur_match = check_how_many_mismatches(	ref_ptr + (left_end_index_in_dir_gen>>1),
																	pattern_ptr + (left_index>>1),					
																	length, true, allowed_mismatches_for_the_current_segment);
				else                //nieparzysta pozycja w tekscie, parzysty offset
					cur_match = check_how_many_mismatches(	ref_ptr + (left_end_index_in_dir_gen>>1),
																	pattern_sft_ptr + (left_index>>1),					
																	length, true, allowed_mismatches_for_the_current_segment);
				break;
			} //end_of_switch

			if((!cur_match) || (cur_match > allowed_mismatches_for_the_current_segment))
				return max_mismatches + 1;
		
			done_mismatches += cur_match;
			allowed_mismatches_for_the_current_segment = allowed_mismatches_for_the_current_segment - cur_match + 1;
			//KOMENTARZ: allowed_mismatches_for_the_current_segment mo¿e byæ unsigned poniewa¿ gdyby cur_match by³o za du¿e
			//to nast¹pi wyjœcie 3 wiersze wy¿ej

			left_offset = CUR_OFFSET - right_index;
			left_index = right_index;
		} //end_for
	}	

	//-------------------------------------------------------------------------------------------------
	//check right part of the read, on the right from a matching part
	//---------------------------------------------------------------------------------------------
	if((done_mismatches <= max_mismatches))
	{
		if(length_after_matching_part)
		{
			uchar_t temp_value = 0;

			temp_value = (sa_part[SA_index] & 1);
			temp_value = (temp_value << 1) | (CUR_OFFSET & 1);
			temp_value = (temp_value << 1) | (predef_match_length & 1);

			switch(temp_value)
			{
			case 0:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index] + predef_match_length) / 2,
																pattern_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1),					
																length_after_matching_part, false, max_mismatches); break;
			case 1:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1),
																pattern_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1),					
																length_after_matching_part + 1, false, max_mismatches); break;
			case 2:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1),
																pattern_sft_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1) + 1,					
																length_after_matching_part, false, max_mismatches); break;
			case 3:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1),
																pattern_sft_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1) + 1,					
																length_after_matching_part + 1, false, max_mismatches); break;
			case 4:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1),
																pattern_sft_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1),					
																length_after_matching_part + 1, false, max_mismatches); break;
			case 5:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1) + 1  ,
																pattern_sft_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1) + 1,					
																length_after_matching_part, false, max_mismatches); break;
			case 6:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1),
																pattern_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1),					
																length_after_matching_part + 1, false, max_mismatches); break;
			case 7:
				cur_match = check_how_many_mismatches(ref_ptr + (sa_part[SA_index]>>1) + (predef_match_length>>1) + 1,
																pattern_ptr + (CUR_OFFSET>>1) + (predef_match_length>>1) + 1,					
																length_after_matching_part, false, max_mismatches); break;
			}  //end_of_switch
			done_mismatches += cur_match;

		} //end_if(length_after_matching_part)	
	}
	
	return done_mismatches;
};

//*********************************************************************************************************
uint32_t CMappingCore::test_selected_rc_data(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len, uint32_t* sa_part, uint32_t SA_index, uint32_t max_mismatches)
{
	uint32_t length_of_the_right_slice = CUR_OFFSET;
	//ma tak zostac dla RZD:
	uint32_t length_of_the_left_slice  = pattern_len - predef_match_length - CUR_OFFSET;
	
	uint32_t cur_match = 0;
	uint32_t done_mismatches = 0;

	uint32_t left_end_index_of_the_right_slice = ref_size - sa_part[SA_index];
	//ma tak zostac dla RZD:
	uint32_t left_end_index_of_the_left_slice  = left_end_index_of_the_right_slice - pattern_len + CUR_OFFSET;
	
	//ma tak zostac dla RZD:
	uint32_t offset_slice_left_index = pattern_len - CUR_OFFSET;
	//offset_slice_left_index - przesuniecie we wzorcu liczone od lewej. wyznacza gdzie sie zaczyna lewy koniec prawego kawalka
	
	//-----------------------------------------------------------------------------------------------
	//-------------------------------RC GENOME ----------------------------------------------------
	//----------------------------------------------------------------------------------------------
	if(length_of_the_right_slice)
	{

		//chyba niepotrzebne:	uint32_t left_index = 0;
		//uint32_t left_offset = CUR_OFFSET - left_index;
		uint32_t left_offset = CUR_OFFSET;

		//one mismatch per segment is allowed for all segments placed on the right of the current one.
		//That gives: (stage_minor - 1)* 1 mismatches
		uint32_t allowed_mismatches_for_the_current_segment = max_mismatches - (stage_minor - 1);

		for(int32_t i = stage_minor - 1; i >= 0; i--)
		{
			uint32_t length = segment_segments[i].second - segment_segments[i].first;
	
			switch (left_end_index_of_the_right_slice & 1)
			{
			case 0:				//je¿eli w tekscie czytamy od parzystej pozycji (czytajac do przodu)
				if(offset_slice_left_index & 1)	//nieparzysty numer dopasowywanego znaku we wzorcu
					cur_match = check_how_many_mismatches(	ref_ptr + (left_end_index_of_the_right_slice>>1),
																		pattern_sft_ptr + (offset_slice_left_index>>1) + 1,					
																		length, false, allowed_mismatches_for_the_current_segment);
				
				else                //parzysty numer dopasowywanego znaku we wzorcu
					cur_match = check_how_many_mismatches(ref_ptr + (left_end_index_of_the_right_slice>>1), 
																	pattern_ptr + (offset_slice_left_index>>1), 
																	length, false, allowed_mismatches_for_the_current_segment);
				break;

			case 1:				//je¿eli w tekscie czytamy od nieparzystej pozycji (czytajac do przodu)
				if(offset_slice_left_index & 1)	//nieparzysty numer dopasowywanego znaku we wzorcu
					cur_match = check_how_many_mismatches(ref_ptr + (left_end_index_of_the_right_slice>>1), 
																	pattern_ptr + (offset_slice_left_index>>1), 
																	length, true, allowed_mismatches_for_the_current_segment);
																	
				
				else                //parzysty numer dopasowywanego znaku we wzorcu
					cur_match = check_how_many_mismatches(ref_ptr + (left_end_index_of_the_right_slice>>1), 
																	pattern_sft_ptr + (offset_slice_left_index>>1), 
																	length, true, allowed_mismatches_for_the_current_segment);
				break;
			} //end of switch;

			if((!cur_match) || (cur_match > allowed_mismatches_for_the_current_segment))
				return max_mismatches + 1;

			done_mismatches += cur_match;
			allowed_mismatches_for_the_current_segment = allowed_mismatches_for_the_current_segment - cur_match + 1;
			left_end_index_of_the_right_slice += length;
			offset_slice_left_index += length;
		} //end_for
	}	//end of if(length_of_the_right_slice)
	
	if((done_mismatches <= max_mismatches))
	{		
		if(length_of_the_left_slice)
		{
			if(left_end_index_of_the_left_slice & 1)
			
				cur_match = check_how_many_mismatches(ref_ptr + left_end_index_of_the_left_slice / 2 ,
																	  pattern_sft_ptr,
																	  length_of_the_left_slice, true, max_mismatches);
			else
				cur_match = check_how_many_mismatches(ref_ptr + left_end_index_of_the_left_slice / 2,
																	  pattern_ptr,	
																	  length_of_the_left_slice, false, max_mismatches);
			done_mismatches += cur_match;
		}
	}       //end_of_if((done_mismatches <= max_mismatches))
	
	return done_mismatches;
};

// ************************************************************************************
uint32_t CMappingCore::check_how_many_mismatches(uchar_t* T, uchar_t* P, uint32_t curr_pattern_len, bool odd_pos, uint32_t max_mismatches)
{
	uint32_t no_of_mismatches		   = 0;
	uint32_t i						   = 0;
	uint32_t curr_pattern_len_in_full_bytes = curr_pattern_len / 2;

	switch (odd_pos)
	{
	case true:
		if(curr_pattern_len & 1)	//nieparzysta pozycja w tekscie i nieparzysta liczba znaków
		{                   //czyli mamy pó³ bajta z przodu
			no_of_mismatches += LO_NIBBLE(T[i] ^ P[i]) != 0;
			if(no_of_mismatches > max_mismatches)
					break;
			i++;
			for(; i <= curr_pattern_len_in_full_bytes; i++ )
			{
				no_of_mismatches += cmp_lut_ptr[T[i] ^ P[i]];
				if(no_of_mismatches > max_mismatches)
					break;
			}
		}
		else		//czyli nieparzysta pozycja w tekscie i parzysta liczba znaków
		{
			no_of_mismatches += LO_NIBBLE(T[i] ^ P[i]) != 0;
			if(no_of_mismatches > max_mismatches)
					break;
			i++;

			for(; i < (curr_pattern_len_in_full_bytes); i++ )
			{
				no_of_mismatches += cmp_lut_ptr[T[i] ^ P[i]];
				if(no_of_mismatches > max_mismatches)
					break;
			}
			if(no_of_mismatches <= max_mismatches)
				no_of_mismatches += HI_NIBBLE(T[i] ^ P[i]) != 0;
		 }
		break;
	case false:
		if(curr_pattern_len & 1)	//parzysta pozycja w tekscie i nieparzysta liczba znaków
		{
			for(i = 0; i < curr_pattern_len_in_full_bytes; i++ )
			{
				no_of_mismatches += cmp_lut_ptr[T[i] ^ P[i]];
				if(no_of_mismatches > max_mismatches)
					break;
			}
			if(no_of_mismatches <= max_mismatches)
				no_of_mismatches += HI_NIBBLE(T[i] ^ P[i]) != 0;
		}
		else		//czyli parzysta pozycja w tekscie i parzysta liczba znaków
		{
			for(i = 0; i < (curr_pattern_len_in_full_bytes); i++ )
			{
				no_of_mismatches += cmp_lut_ptr[T[i] ^ P[i]];

				if(no_of_mismatches > max_mismatches)
					break;
			}
		}
		break;
	}
	return no_of_mismatches;
};

uint32_t CMappingCore::calculate_CF_for_dir_and_rc_gen(uint32_t pattern_len, uint32_t position_read_from_SA_array, uint32_t fixedLeft, uint32_t fixedRight, bool dir_gen_flag)
{
	uint32_t left_end;
	if(dir_gen_flag)
		left_end = position_read_from_SA_array - CUR_OFFSET;
	else
		left_end = ref_size - position_read_from_SA_array - pattern_len + CUR_OFFSET;
	uint32_t CF_value = 0;
	uchar_t* T_ptr = ref_ptr + (left_end >> 1);
	uint32_t i;
	uint32_t len_in_bytes;

	//---------------------CF dla czesci przed matching part:
	if(fixedLeft)
	{
		fixedLeft++;
		len_in_bytes = fixedLeft >> 1;

		switch(left_end & 1)
		{
		case 0:	//parzysta pozycja w tekscie
			for(i = 0; i < len_in_bytes; i++)
				CF_value += cf_lut[*(T_ptr + i)];

			//if(fixedLeft & 1)		//je¿eli liczba znaków w czesci przed "matching part" jest nieparzysta
				//CF_value += cf_lut[RAW_HI_NIBBLE(*(T_ptr + i))] - 0x01000000ul ;
			break;
		case 1: //nieparzysta pozycja w tekscie
			i = 0;
			CF_value += cf_lut[LO_NIBBLE(*(T_ptr + i))] - 0x01000000ul;
		
			for(i = 1; i < len_in_bytes; i++)
				CF_value += cf_lut[*(T_ptr + i)];

			//if(fixedLeft & 1)		//je¿eli liczba znaków w czesci przed "matching part" jest nieparzysta
				//CF_value += cf_lut[*(T_ptr + i)];
			//else
				CF_value += cf_lut[RAW_HI_NIBBLE(*(T_ptr + i))] - 0x01000000ul ;
		}// end_of_switch
	}

	if(fixedRight != pattern_len)
	{
		if(fixedRight & 1)
			fixedRight--;			//fixedRight is always even;
		
		uint32_t len_in_chars = pattern_len - fixedRight;
		len_in_bytes = len_in_chars >> 1;

		left_end += fixedRight;
		T_ptr = ref_ptr + (left_end >> 1);

		switch(left_end & 1)
		{
		case 0:	//parzysta pozycja w tekscie
			for(i = 0; i < len_in_bytes; i++)
			{
				CF_value += cf_lut[*(T_ptr + i)];
				
			}
			if((len_in_chars & 1))		//je¿eli liczba znaków w readzie jest nieparzysta
				//CF_value += cf_lut[*(T_ptr + i)];
				CF_value += cf_lut[RAW_HI_NIBBLE(*(T_ptr + i))] - 0x01000000ul ;
			break;
		case 1: //nieparzysta pozycja w tekscie
			i = 0;
			CF_value += cf_lut[LO_NIBBLE(*(T_ptr + i))] - 0x01000000ul;
		
			for(i = 1; i < len_in_bytes; i++)
				CF_value += cf_lut[*(T_ptr + i)];

			if(len_in_chars & 1)		//je¿eli liczba znaków w readzie jest nieparzysta
				CF_value += cf_lut[*(T_ptr + i)];
			else
				CF_value += cf_lut[RAW_HI_NIBBLE(*(T_ptr + i))] - 0x01000000ul ;
		} //end_of_switch

	}
		

	return CF_value;
}

//*********************************************************************************************
uint32_t CMappingCore::calculate_CF_for_a_pattern(uchar_t *pattern, uint32_t pattern_len, uint32_t fixedLeft, uint32_t fixedRight)
{
	uint32_t CF_value = 0;
	uint32_t i, j;
	uint32_t len_in_bytes;

	//---------------------CF dla czesci przed matching part:
	if(fixedLeft)
	{
		fixedLeft++;
		len_in_bytes = fixedLeft >> 1;

		for(i = 0; i < len_in_bytes; i++)
			CF_value += cf_lut[ pattern[i] ];
	}
	//--------------------policzony kawa³ak przed matching part
	//--------------------teraz bedzie "za" matching part

	if(fixedRight != pattern_len)
	{
		if(fixedRight & 1)
			fixedRight--;			//fixedRight is always even

	
		uint32_t len_in_chars = pattern_len - fixedRight;
		len_in_bytes = len_in_chars >> 1;

		for(i = (fixedRight >> 1), j = 0; j < len_in_bytes; i++, j++)
					CF_value += cf_lut[ pattern[i] ];


		if(len_in_chars & 1)		//je¿eli liczba znaków w czesci za "matching part" jest nieparzysta
			CF_value += cf_lut[ RAW_HI_NIBBLE(pattern[i]) ] - 0x01000000ul ;
	}

	return CF_value;
}

//*********************************************************************************************
#define ABS_DIFF(x, y)	((x) > (y) ? (x) - (y) : (y) - (x))

bool CMappingCore::test_CF_difference(uint32_t CF_1, uint32_t CF_2, uint32_t max_diff)
{
	uint32_t CF_diff;
	uint32_t a, b;

	a = CF_1 & 0xff;
	b = CF_2 & 0xff;
	CF_diff = ABS_DIFF(a, b);
	if(CF_diff > max_diff)
		return false;
	CF_1 >>= 8;
	CF_2 >>= 8;

	a = CF_1 & 0xff;
	b = CF_2 & 0xff;
	CF_diff += ABS_DIFF(a, b);
	if(CF_diff > max_diff)
		return false;
	CF_1 >>= 8;
	CF_2 >>= 8;

	a = CF_1 & 0xff;
	b = CF_2 & 0xff;
	CF_diff += ABS_DIFF(a, b);
	if(CF_diff > max_diff)
		return false;
	CF_1 >>= 8;
	CF_2 >>= 8;

	CF_diff += ABS_DIFF(CF_1, CF_2);

	return CF_diff <= max_diff;
}

#endif

// EOF
