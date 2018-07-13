// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.1
// Date    : 2018-07-10
// License : GNU GPL 3
// *******************************************************************************************


#include "mapping_core.h"
#include <algorithm>
#include <utility>
#include <cmath>
#include "../libs/asmlib.h"
#include "../common/utils.h"
#include <assert.h>

//**********************************************************************************************************
bool  CMappingCore::findExactMatch_diff_len_version(uchar_t* pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
	bool basic_data_part_is_new,
	uint32_t *sa_last_index, uint32_t* sa_start_index,
	uint32_t *sa_part, uint32_t sa_size,
	bool dir_gen_flag)
{
	int cmp = 1;
	uint32_t cur_SA_index;

	if (basic_data_part_is_new)
	{
		cur_SA_index = *sa_last_index;
		if (cur_SA_index < sa_size)
		{
			do
			{
				//checking identity at length "min_read_length" 
				//"pattern_len" parameter is needed for correct calculations for rc
				cmp = check_basic_part_of_diff_len_reads(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part[cur_SA_index], dir_gen_flag);

				cur_SA_index++;
			} while ((cmp < 0) && (cur_SA_index < sa_size));

			cur_SA_index--;			//in this way (cur_SA_index) contains the last checked position
			if (cmp == 0)
			{
				*sa_start_index = cur_SA_index;

				cmp = search_in_unknown_range_for_exact_match_diff_len_version(pattern_ptr, pattern_sft_ptr, pattern_len,
					sa_part, sa_size,
					*sa_start_index, sa_last_index, dir_gen_flag, 0);
			}
			else
			{
				*sa_start_index = cur_SA_index;
				*sa_last_index = cur_SA_index;
			}
		}
		else //if(cur_SA_index < sa_size)
		{
			*sa_start_index = cur_SA_index;
		}
	}
	else
		cmp = search_in_known_range_for_exact_match_diff_len_version(pattern_ptr, pattern_sft_ptr, pattern_len,
			sa_part, sa_size,
			*sa_start_index, *sa_last_index, dir_gen_flag, 0);

	return (cmp == 0);	//exact_match_found; 
}

// ************************************************************************************
uint32_t CMappingCore::findMismatches_dir_and_rc(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
	bool matching_part_is_new,
	uint32_t *sa_part, uint32_t sa_size,
	uint32_t* sa_last_index, uint32_t* sa_start_index,
	bool dir_gen_flag, bool* inside_malicious_group, bool* aux_struct_valid, uint32_t max_mismatches)
{
	uint32_t cur_SA_index;
	uint32_t no_of_mismatches = max_mismatches + 1;
	int32_t  cmp;

	uint32_t cur_malicious_group_length = malicious_group_length;
	if (sensitive_mode)
		cur_malicious_group_length *= sens_malicious_group_factor;

	if (matching_part_is_new) //"matching_part" appeared for the first time
	{
		cur_SA_index = *sa_last_index;
		if (cur_SA_index < sa_size)
		{
#ifdef NO_BINARY_SEARCH
			do
			{
				cmp = check_matching_part(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part[cur_SA_index], dir_gen_flag);
				cur_SA_index++;

			} while ((cmp < 0) && (cur_SA_index < sa_size));

			cur_SA_index--;
#else
			cmp = check_matching_part(pattern_ptr, pattern_sft_ptr, pattern_len, sa_part[cur_SA_index], dir_gen_flag);
			if (cmp < 0)
			{
				uint32_t move_size = 0;
				cmp = binary_sa_search(pattern_ptr, pattern_sft_ptr, pattern_len, predef_match_length, sa_part + cur_SA_index, sa_size - cur_SA_index, &move_size, dir_gen_flag);
				cur_SA_index += move_size;
			}
#endif

			if (cmp == 0)
			{
				*sa_start_index = cur_SA_index;
#ifdef USE_128b_SSE_CF_AND_DP
				if (dir_gen_flag)
					calculate_CF_128b_SSE_for_a_pattern(pattern_ptr, pattern_sft_ptr, min_read_len, start_pos, end_pos, CF_value_SSE_for_a_dir_pattern);
				else
				{
					uint32_t diff_len = pattern_len - min_read_len;

					if (diff_len & 1)	                    //diff_len is odd
						calculate_CF_128b_SSE_for_a_pattern(pattern_sft_ptr + (diff_len >> 1) + 1, pattern_ptr + (diff_len >> 1), min_read_len, min_read_len - end_pos, min_read_len - start_pos, CF_value_SSE_for_a_rc_pattern);
					else                                    //diff_len is even
						calculate_CF_128b_SSE_for_a_pattern(pattern_ptr + (diff_len >> 1), pattern_sft_ptr + (diff_len >> 1), min_read_len, min_read_len - end_pos, min_read_len - start_pos, CF_value_SSE_for_a_rc_pattern);
				}
#endif					
				no_of_mismatches = search_in_unknown_range_dir_and_rc(pattern_ptr, pattern_sft_ptr, pattern_len,
					sa_part, sa_size,
					*sa_start_index, sa_last_index, dir_gen_flag, max_mismatches);
				if ((*sa_last_index - *sa_start_index) >= cur_malicious_group_length)
				{
					*inside_malicious_group = true;
					*aux_struct_valid = false;
				}
				else
					*inside_malicious_group = false;
			}
			else
			{
				*sa_start_index = cur_SA_index;
				*sa_last_index = cur_SA_index;

				*inside_malicious_group = false;
			}

		} //end_if(cur_SA_index < sa_size)
		  //else
		  //the end of the array, the pattern was not found
	}
	else		//matching_part_is_new == false
	{

#ifdef USE_128b_SSE_CF_AND_DP
		if (dir_gen_flag)
			calculate_CF_128b_SSE_for_a_pattern(pattern_ptr, pattern_sft_ptr, min_read_len, start_pos, end_pos, CF_value_SSE_for_a_dir_pattern);
		else
		{
			uint32_t diff_len = pattern_len - min_read_len;

			if (diff_len & 1)						//diff_len is odd
				calculate_CF_128b_SSE_for_a_pattern(pattern_sft_ptr + (diff_len >> 1) + 1, pattern_ptr + (diff_len >> 1), min_read_len, min_read_len - end_pos, min_read_len - start_pos, CF_value_SSE_for_a_rc_pattern);
			else                                   //diff_len is even
				calculate_CF_128b_SSE_for_a_pattern(pattern_ptr + (diff_len >> 1), pattern_sft_ptr + (diff_len >> 1), min_read_len, min_read_len - end_pos, min_read_len - start_pos, CF_value_SSE_for_a_rc_pattern);
		}
#endif
		if (*inside_malicious_group)
			if (dir_gen_flag)
#ifdef USE_128b_SSE_CF_AND_DP
				no_of_mismatches = search_in_malicious_group_with_DP(pattern_ptr, pattern_sft_ptr, pattern_len,
					sa_part + *sa_start_index, *sa_last_index - *sa_start_index,
					dir_aux_arr, dir_aux_params, &dir_aux_struct_valid, true, max_mismatches,
					*sa_start_index);

#else
				no_of_mismatches = search_in_malicious_group_dir_and_rc(pattern_ptr, pattern_sft_ptr, pattern_len,
					sa_part + *sa_start_index, *sa_last_index - *sa_start_index,
					dir_aux_arr, dir_aux_params, &dir_aux_struct_valid, true, max_mismatches,
					*sa_start_index);
#endif
			else
				//rc_gen
#ifdef USE_128b_SSE_CF_AND_DP
				no_of_mismatches = search_in_malicious_group_with_DP(pattern_ptr, pattern_sft_ptr, pattern_len,
					sa_part + *sa_start_index, *sa_last_index - *sa_start_index,
					rc_aux_arr, rc_aux_params, &rc_aux_struct_valid, false, max_mismatches,
					*sa_start_index);
#else
				no_of_mismatches = search_in_malicious_group_dir_and_rc(pattern_ptr, pattern_sft_ptr, pattern_len,
					sa_part + *sa_start_index, *sa_last_index - *sa_start_index,
					rc_aux_arr, rc_aux_params, &rc_aux_struct_valid, false, max_mismatches,
					*sa_start_index);
#endif
		else

			no_of_mismatches = search_in_known_range_dir_and_rc(pattern_ptr, pattern_sft_ptr, pattern_len,
				sa_part, sa_size,
				*sa_start_index, *sa_last_index, dir_gen_flag, max_mismatches);

	}

	return no_of_mismatches;
}

// EOF 