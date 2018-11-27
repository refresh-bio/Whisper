// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : see defs.h
// Date    : see defs.h
// License : GNU GPL 3
// *******************************************************************************************


#include "mapping_core.h"
#include <algorithm>
#include <utility>
#include <assert.h>

//#define LEV_MYERS_ASSERTIONS

//**********************************************************************************************************
uint32_t CMappingCore::search_in_unknown_range_for_exact_match_diff_len_version(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
	uint32_t* sa_part, uint32_t sa_size, uint32_t start_SA_index, uint32_t* last_index,
	bool dir_gen_flag, uint32_t max_mismatches)
	//pattern_len is only passed on as a parameter
{
	uint32_t i = start_SA_index;
	uint32_t cur_match = 0;
	uint32_t best_match = max_mismatches + 1;
	int32_t  cmp = 0;

	if (pattern_len == min_read_len)
	{
		do
		{
			//There is nothing more to check. The basic_part contains the whole read and it has already been checked
			best_match = 0;
			if (dir_gen_flag)
				dir_match_positions.push_back(make_pair(sa_part[i], cur_match));
			else
				rc_match_positions.push_back(make_pair(ref_size - (sa_part[i]) - 1, cur_match));

			i++;
			if (i < sa_size)
				cmp = check_basic_part_of_diff_len_reads(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part[i], dir_gen_flag);
		} while ((i < sa_size) && (cmp <= 0));
	}
	else
	{
		// The read is longer than basic_part and its ending must be checked
		while ((i < sa_size) && (cmp <= 0))
		{
			//-------------------------------DIR GENOME ----------------------------------------
			if (dir_gen_flag)
			{
				auto even_min_read_len = min_read_len & ~1u;

				switch (sa_part[i] & 1)
				{
				case 1:
					cur_match = compare_dir_str_with_gen(ref_ptr + (sa_part[i] + even_min_read_len) / 2, pattern_sft_ptr + even_min_read_len / 2, pattern_len - even_min_read_len, true);
					break;
				case 0:
					cur_match = compare_dir_str_with_gen(ref_ptr + (sa_part[i] + even_min_read_len) / 2, pattern_ptr + even_min_read_len / 2, pattern_len - even_min_read_len, false);
					break;
				}
			}
			//--------------------------------RC GENOME ---------------------------------------
			else
			{
				uint32_t right_end_index_in_dir_gen;				                    //starting point (when reading backward) in the reference 	
				right_end_index_in_dir_gen = ref_size - sa_part[i] - 1 - min_read_len;

				uint32_t matching_part_right_point = pattern_len - 1 - min_read_len;	//starting point (when reading backward) in the pattern

				switch (right_end_index_in_dir_gen & 1)
				{
				case 0:				                //we read from an even position in the text (when reading backward)
					if (matching_part_right_point & 1)	//The odd number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1) + 1, pattern_len - min_read_len, false);

					else                                 //Even number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), pattern_len - min_read_len, false);
					break;

				case 1:				                 //we read from an odd position in the text (when reading backward)
					if (matching_part_right_point & 1)	 //The odd number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), pattern_len - min_read_len, true);
					else                                 //Even number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1), pattern_len - min_read_len, true);
					break;
				} //end of switch;
			} //end RC_GENOME
			  //-------------------------------------------------------------------------------

			if (cur_match == 0)
			{
				best_match = 0;
				if (dir_gen_flag)
				{
					dir_match_positions.push_back(make_pair(sa_part[i], cur_match));
				}
				else
				{
					rc_match_positions.push_back(make_pair(ref_size - (sa_part[i]) - 1, cur_match));
				}
			}

			i++;
			if (i < sa_size)
			{
				cmp = check_basic_part_of_diff_len_reads(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part[i], dir_gen_flag);
			}
		} //end_while(( i < sa_size) && (cmp <= 0))
	}
	*last_index = i;

	return best_match;
}

// ************************************************************************************ 22.10.2017 ok
uint32_t CMappingCore::search_in_known_range_for_exact_match_diff_len_version(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
	uint32_t* sa_part, uint32_t sa_size, uint32_t start_SA_index, uint32_t last_index,
	bool dir_gen_flag, uint32_t max_mismatches)
{
	uint32_t i = start_SA_index;
	uint32_t cur_match = 0;
	uint32_t best_match = max_mismatches + 1;

	// TODO (Future): Here we can probably assume that the basic_part is already matched and check only the endings
	if (start_SA_index < last_index)
	{
		while ((i < last_index) && (i < sa_size))
		{
			//-------------------------------DIR GENOME ----------------------------------------------------
			if (dir_gen_flag)
			{
				switch (sa_part[i] & 1)
				{
				case 1:				//we read from an odd position in the text
					cur_match = compare_dir_str_with_gen(ref_ptr + sa_part[i] / 2, pattern_sft_ptr, pattern_len, true);
					break;

				case 0:				//we read from an even position in the text
					cur_match = compare_dir_str_with_gen(ref_ptr + sa_part[i] / 2, pattern_ptr, pattern_len, false);
					break;
				} //end_of_switch
			}
			//--------------------------------RC GENOME -----------------------------------------------------
			else
			{
				uint32_t right_end_index_in_dir_gen;				    //starting point (when reading backward) in the reference 	
				right_end_index_in_dir_gen = ref_size - sa_part[i] - 1;

				uint32_t matching_part_right_point = pattern_len - 1;	//starting point (when reading backward) in the pattern

				switch (right_end_index_in_dir_gen & 1)
				{
				case 0:				                //we read from an even position in the text (when reading backward)
					if (matching_part_right_point & 1)	//The odd number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1) + 1, pattern_len, false);

					else                                //Even number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), pattern_len, false);
					break;

				case 1:				                //we read from an odd position in the text (when reading backward)
					if (matching_part_right_point & 1)	//The odd number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), pattern_len, true);

					else                                //Even number of the character being matched in the pattern
						cur_match = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1), pattern_len, true);
					break;
				} //end of switch;
			}
			//-----------------------------------------------------------------------------------------------

			if (cur_match == 0)
			{
				best_match = 0;
				if (dir_gen_flag)
				{
					dir_match_positions.push_back(make_pair(sa_part[i], cur_match));
				}
				else
				{
					rc_match_positions.push_back(make_pair(ref_size - (sa_part[i]) - 1, cur_match));
				}
			}
			i++;
		}
		return best_match;
	}

	return max_mismatches + 1;
}

//**********************************************************************************************************
//Traversing a suffix array when trying to match a pattern. 
//The traversing starts at start_SA_index, last_index is unknown (output parameter)
//***********************************************************************************************************
uint32_t CMappingCore::search_in_unknown_range_dir_and_rc(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
	uint32_t* sa_part, uint32_t sa_size, uint32_t start_SA_index, uint32_t* last_index,
	bool dir_gen_flag, uint32_t max_mismatches)
	//pattern_len is only passed on as a parameter
{
	uint32_t i = start_SA_index;
	uint32_t cur_match = 0;
	uint32_t best_match = max_mismatches + 1;
	int32_t  cmp = 0;

	LevMyers* levMyers;
	levMyers = (pattern_len < 256)
		? (pattern_len < 128) ? levMyers128 : levMyers256
		: levMyers64;

	if (max_mismatches >= lev_alg_thr)
	{
#ifdef LEV_MYERS_ASSERTIONS
		edit_dist->LevMyers_PP(pattern_ptr, pattern_len, genome_t::direct);		// Direction is always direct as read is already reverse-complemented (if dir_falg is false)
#endif
		levMyers->preprocess(pattern_ptr, pattern_len, genome_t::direct);
	}
	while ((i < sa_size) && (cmp <= 0))
	{
		//-------------------------------DIR GENOME ----------------------------------------
		if (dir_gen_flag)
		{
#ifdef USE_128b_SSE_CF_AND_DP
			calculate_CF_128b_SEE_for_dir_and_rc_gen(start_pos, end_pos, min_read_len, sa_part[i], CF_vector_128b_SSE_dir_gen[i], dir_gen_flag);

			if (test_CF_128b_SSE_difference(CF_value_SSE_for_a_dir_pattern, CF_vector_128b_SSE_dir_gen[i], 4 * max_mismatches))
			{
				no_of_CF_accepted++;
				ref_pos_t ref_pos;
				if (max_mismatches >= lev_alg_thr) {
					levMyers->dynamicProgramming(sa_part[i] - start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos, cur_match);
#ifdef LEV_MYERS_ASSERTIONS
					uint32_t ref_pos2, cur_match2;
					edit_dist->LevMyers(sa_part[i] - start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos2, cur_match2);
					ASSERT(ref_pos == ref_pos2 && cur_match == cur_match2, "CMappingCore::search_in_unknown_range_dir_and_rc() error!");
#endif 
				}
				else
					cur_match = edit_dist->LevDiag(pattern_ptr, pattern_len, sa_part[i], start_pos, end_pos, max_mismatches, 1);

				if (cur_match <= max_mismatches)
					no_of_Lev_positive++;
			}
			else
			{
				cur_match = max_mismatches + 1;
				no_of_CF_discarded++;
			}

#else
			cur_match = test_selected_data(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part, i, max_mismatches);
#endif
		}

		//--------------------------------RC GENOME ---------------------------------------
		else
		{
#ifdef USE_128b_SSE_CF_AND_DP
			calculate_CF_128b_SEE_for_dir_and_rc_gen(min_read_len - end_pos, min_read_len - start_pos, min_read_len, sa_part[i], CF_vector_128b_SSE_rc_gen[i], 0);

			if (test_CF_128b_SSE_difference(CF_value_SSE_for_a_rc_pattern, CF_vector_128b_SSE_rc_gen[i], 4 * max_mismatches))
			{
				ref_pos_t ref_pos;
				if (max_mismatches >= lev_alg_thr)
				{
					levMyers->dynamicProgramming(ref_size - sa_part[i] - pattern_len + start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos, cur_match);
#ifdef LEV_MYERS_ASSERTIONS
					uint32_t ref_pos2, cur_match2;
					edit_dist->LevMyers(ref_size - sa_part[i] - pattern_len + start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos2, cur_match2);
					ASSERT(ref_pos == ref_pos2 && cur_match == cur_match2, "CMappingCore::search_in_unknown_range_dir_and_rc() error!");
#endif
				}
				else
					cur_match = edit_dist->LevDiag(pattern_ptr, pattern_len, sa_part[i], pattern_len - end_pos, pattern_len - start_pos, max_mismatches, 0);

				if (cur_match <= max_mismatches)
					no_of_Lev_positive++;
			}
			else
			{
				cur_match = max_mismatches + 1;
				no_of_CF_discarded++;
			}
#else
			cur_match = test_selected_rc_data(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part, i, max_mismatches);
#endif
		}	//end_of RC_Genome
			//--------------------------------------------------------------------------------

		if (cur_match < best_match)
			best_match = cur_match;

		// Future: consider saving all found mappings not just best best one (important for all strata and partially also for first stratum)
		if (cur_match <= max_mismatches)
		{
			if (dir_gen_flag)
			{
				dir_match_positions.push_back(make_pair(sa_part[i] - CUR_OFFSET, cur_match));
			}
			else
			{
				rc_match_positions.push_back(make_pair(ref_size - (sa_part[i] - CUR_OFFSET) - 1, cur_match));
			}
		}

		i++;
		if (i < sa_size)
		{
			cmp = check_matching_part(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part[i], dir_gen_flag);
		}
	}
	*last_index = i;

	return best_match;
}

//**********************************************************************************************************
//Traversing a suffix array when trying to match a pattern.
//The traversing starts at start_SA_index, stops at last_index 
//**********************************************************************************************************
uint32_t CMappingCore::search_in_known_range_dir_and_rc(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
	uint32_t* sa_part, uint32_t sa_size, uint32_t start_SA_index, uint32_t last_index,
	bool dir_gen_flag, uint32_t max_mismatches)
{
	uint32_t i = start_SA_index;
	uint32_t cur_match = 0;
	uint32_t best_match = max_mismatches + 1;

	LevMyers* levMyers;
	levMyers = (pattern_len < 256)
		? (pattern_len < 128) ? levMyers128 : levMyers256
		: levMyers64;

	if (max_mismatches >= lev_alg_thr)
	{
#ifdef LEV_MYERS_ASSERTIONS
		edit_dist->LevMyers_PP(pattern_ptr, pattern_len, genome_t::direct);		// Direction is always direct as read is already reverse-complemented (if dir_falg is false)
#endif
		levMyers->preprocess(pattern_ptr, pattern_len, genome_t::direct);
	}
	if (start_SA_index < last_index)
	{
		while ((i < last_index) && (i < sa_size))
		{
			//-------------------------------DIR GENOME ----------------------------------------------------
			if (dir_gen_flag)
			{
#ifdef USE_128b_SSE_CF_AND_DP
				if (test_CF_128b_SSE_difference(CF_value_SSE_for_a_dir_pattern, CF_vector_128b_SSE_dir_gen[i], 4 * max_mismatches))
				{
					no_of_CF_accepted++;
					ref_pos_t ref_pos;
					if (max_mismatches >= lev_alg_thr)
					{
						levMyers->dynamicProgramming(sa_part[i] - start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos, cur_match);
#ifdef LEV_MYERS_ASSERTIONS
						uint32_t ref_pos2, cur_match2;
						edit_dist->LevMyers(sa_part[i] - start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos2, cur_match2);
						ASSERT(ref_pos == ref_pos2 && cur_match == cur_match2, "CMappingCore::search_in_known_range_dir_and_rc() error!");
#endif
					}
					else
						cur_match = edit_dist->LevDiag(pattern_ptr, pattern_len, sa_part[i], start_pos, end_pos, max_mismatches, 1);

					if (cur_match <= max_mismatches)
						no_of_Lev_positive++;
				}
				else
				{
					no_of_CF_discarded++;
					cur_match = max_mismatches + 1;
				}
#else
				cur_match = test_selected_data(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part, i, max_mismatches);
#endif
			}
			//--------------------------------RC GENOME -----------------------------------------------------
			else
			{
#ifdef USE_128b_SSE_CF_AND_DP
				if (test_CF_128b_SSE_difference(CF_value_SSE_for_a_rc_pattern, CF_vector_128b_SSE_rc_gen[i], 4 * max_mismatches))
				{
					ref_pos_t ref_pos;
					if (max_mismatches >= lev_alg_thr) {
						levMyers->dynamicProgramming(ref_size - sa_part[i] - pattern_len + start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos, cur_match);
#ifdef LEV_MYERS_ASSERTIONS
						uint32_t ref_pos2, cur_match2;
						edit_dist->LevMyers(ref_size - sa_part[i] - pattern_len + start_pos - max_mismatches, pattern_len + 2 * max_mismatches, max_mismatches, ref_pos2, cur_match2);
						ASSERT(ref_pos == ref_pos2 && cur_match == cur_match2, "CMappingCore::search_in_known_range_dir_and_rc() error!");
#endif
					}
					else
						cur_match = edit_dist->LevDiag(pattern_ptr, pattern_len, sa_part[i], pattern_len - end_pos, pattern_len - start_pos, max_mismatches, 0);

					if (cur_match <= max_mismatches)
						no_of_Lev_positive++;
				}
				else
				{
					no_of_CF_discarded++;
					cur_match = max_mismatches + 1;
				}
#else
				cur_match = test_selected_rc_data(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part, i, max_mismatches);
#endif
			}//end_of RC_Genome
			 //-----------------------------------------------------------------------------------------------

			if (cur_match < best_match)
				best_match = cur_match;
			if (cur_match <= max_mismatches)
			{
				if (dir_gen_flag)
				{
					dir_match_positions.push_back(make_pair(sa_part[i] - CUR_OFFSET, cur_match));
				}
				else
				{
					rc_match_positions.push_back(make_pair(ref_size - (sa_part[i] - CUR_OFFSET) - 1, cur_match));
				}
			}
			i++;
		}

		return best_match;
	}

	return max_mismatches + 1;
}

// ************************************************************************************
int32_t CMappingCore::binary_sa_search(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len, uint32_t len_to_compare,
	uint32_t* SA, uint32_t sa_size,
	uint32_t *idx, bool dir_gen_flag)
{
	uint32_t size, half;
	uint32_t i;
	int32_t r = -1;
	bool matching_part_found = false;

	uint32_t skip;

	if (dir_gen_flag)
		skip = sa_dir_skip;
	else
		skip = sa_rc_skip;
	i = skip;

	//-----------------------------------------------------------------------------------
	//start skipping phase
	//-----------------------------------------------------------------------------------
	while (i < sa_size)
	{
		r = check_matching_part(pattern_ptr, pattern_sft_ptr, pattern_len, SA[i], dir_gen_flag);

		if (r < 0)
			i += skip;
		else //if(r >= 0) 
		{
			if (r == 0)
				matching_part_found = true;
			break;
		}
	}//end_while;
	 //-----------------------------------------------------------------------------------
	 //end of skipping phase
	 //---------------------------------------------------------------------------------

	if (r < 0)
		size = sa_size + skip - i;
	else
		size = skip;

	//-----------------------------------------------------------------------------------
	//start binary search
	//---------------------------------------------------------------------------------
	for (i = i - skip, half = size >> 1; 0 < size; size = half, half >>= 1)
	{
		r = check_matching_part(pattern_ptr, pattern_sft_ptr, pattern_len, SA[i + half], dir_gen_flag);

		if (r < 0)
		{
			i += half + 1;
			half -= (size & 1) ^ 1;
		}
		else
			if (r == 0)
				matching_part_found = true;
	} //end_for

	*idx = i;
	if (matching_part_found)
		return 0;
	else
		return 1;
}

// EOF 