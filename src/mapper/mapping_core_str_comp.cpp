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
#include "mapping_core_luts.h"

#include <assert.h>

//**********************************************************************************************************
void CMappingCore::string_copy(uchar_t *target, uchar_t *source, uint32_t len)
//len - a number of CHARACTERS in the string
{
	uint32_t no_of_full_bytes = len / 2;
	uint32_t i;

	for (i = 0; i < no_of_full_bytes; i++)
		target[i] = source[i];

	if (len & 1)
		target[i] = (source[i] & 0xf0);
}

//**********************************************************************************************************
int32_t CMappingCore::compare_dir_str_with_gen(uchar_t* text_ptr, uchar_t *pattern_ptr, uint32_t len_to_compare, bool odd_pos)
{
	uint32_t i = 0;
	uint32_t len_to_compare_in_full_bytes = len_to_compare / 2;
	switch (odd_pos)
	{
	case true:
		if (len_to_compare & 1)	//An odd position in the text and an odd number of characters
		{
			if (LO_NIBBLE(text_ptr[i]) < LO_NIBBLE(pattern_ptr[i])) return -1;
			if (LO_NIBBLE(text_ptr[i]) > LO_NIBBLE(pattern_ptr[i])) return 1;
			i++;

			return MEMCMP(text_ptr + 1, pattern_ptr + 1, len_to_compare_in_full_bytes);
		}
		else		//An odd position in the text and an even number of characters
		{
			if (LO_NIBBLE(text_ptr[i]) < LO_NIBBLE(pattern_ptr[i])) return -1;
			if (LO_NIBBLE(text_ptr[i]) > LO_NIBBLE(pattern_ptr[i])) return 1;
			i++;

			int res = MEMCMP(text_ptr + 1, pattern_ptr + 1, len_to_compare_in_full_bytes - 1);
			if (res != 0)
				return res;

			if (HI_NIBBLE(text_ptr[len_to_compare_in_full_bytes]) < HI_NIBBLE(pattern_ptr[len_to_compare_in_full_bytes])) return -1;
			if (HI_NIBBLE(text_ptr[len_to_compare_in_full_bytes]) > HI_NIBBLE(pattern_ptr[len_to_compare_in_full_bytes])) return 1;
		}
		break;
	case false:
		if (len_to_compare & 1)	//even position in the text and an odd number of characters
		{

			int res = MEMCMP(text_ptr, pattern_ptr, len_to_compare_in_full_bytes);
			if (res != 0)
				return res;

			if (HI_NIBBLE(text_ptr[len_to_compare_in_full_bytes]) < HI_NIBBLE(pattern_ptr[len_to_compare_in_full_bytes])) return -1;
			if (HI_NIBBLE(text_ptr[len_to_compare_in_full_bytes]) > HI_NIBBLE(pattern_ptr[len_to_compare_in_full_bytes])) return 1;
		}
		else		//Even position in the text and even number of characters
			return MEMCMP(text_ptr, pattern_ptr, len_to_compare_in_full_bytes);
		break;
	}

	return 0;	//the whole string has been analyzed, which means there is a match
}

//**********************************************************************************************************
int32_t CMappingCore::compare_rc_str_backward(uchar_t* text_ptr, uchar_t *pattern_ptr, uint32_t len_to_compare, bool odd_pos)
{
	int32_t bytes_to_compare = len_to_compare / 2;
	int32_t i;

	switch (odd_pos)
	{
	case true:
		if (len_to_compare & 1)	//An odd position in the text and an odd number of characters
		{
			for (i = 0; i < bytes_to_compare; i++)
			{
				if (*text_ptr != *pattern_ptr)
					return CMP_DIFFERENT_SWAPPED_NIBBLES(*text_ptr, *pattern_ptr);
				text_ptr--; pattern_ptr--;
			}
			if (LO_NIBBLE(*text_ptr) < LO_NIBBLE(*pattern_ptr)) return 1;
			if (LO_NIBBLE(*text_ptr) > LO_NIBBLE(*pattern_ptr)) return -1;
		}
		else		//An odd position in the text and an even number of characters
		{
			for (i = 0; i < bytes_to_compare; i++)
			{
				if (*text_ptr != *pattern_ptr)
					return CMP_DIFFERENT_SWAPPED_NIBBLES(*text_ptr, *pattern_ptr);
				text_ptr--; pattern_ptr--;
			}
		}
		break;
	case false:
		if (len_to_compare & 1)	//even position in the text and an odd number of characters
		{
			if (HI_NIBBLE(*text_ptr) < HI_NIBBLE(*pattern_ptr)) return 1;
			if (HI_NIBBLE(*text_ptr) > HI_NIBBLE(*pattern_ptr)) return -1;
			text_ptr--; pattern_ptr--;

			for (i = 0; i < bytes_to_compare; i++)
			{
				if (*text_ptr != *pattern_ptr)
					return CMP_DIFFERENT_SWAPPED_NIBBLES(*text_ptr, *pattern_ptr);
				text_ptr--; pattern_ptr--;
			}
		}
		else		//Even position in the text and even number of characters
		{
			if (HI_NIBBLE(*text_ptr) < HI_NIBBLE(*pattern_ptr)) return 1;
			if (HI_NIBBLE(*text_ptr) > HI_NIBBLE(*pattern_ptr)) return -1;
			text_ptr--; pattern_ptr--;

			for (i = 1; i < bytes_to_compare; i++)
			{
				if (*text_ptr != *pattern_ptr)
					return CMP_DIFFERENT_SWAPPED_NIBBLES(*text_ptr, *pattern_ptr);
				text_ptr--; pattern_ptr--;;
			}

			if (LO_NIBBLE(*text_ptr) < LO_NIBBLE(*pattern_ptr)) return 1;
			if (LO_NIBBLE(*text_ptr) > LO_NIBBLE(*pattern_ptr)) return -1;
		}
		break;
	}

	return 0;	//the whole string has been analyzed, which means there is a match
}

//**********************************************************************************************************
// Check if two strings are the same
bool CMappingCore::compare_str_to_str(uchar_t *pattern_ptr, uint32_t pattern_len, uchar_t *old_pattern_ptr)
{
	uint32_t pattern_len_in_bytes = pattern_len / 2;

	if (old_pattern_ptr)
	{
		if (MEMCMP(pattern_ptr, old_pattern_ptr, pattern_len_in_bytes))
			return false;
		if (pattern_len & 1)
			return HI_NIBBLE(pattern_ptr[pattern_len_in_bytes]) == HI_NIBBLE(old_pattern_ptr[pattern_len_in_bytes]);
		return true;
	}
	else
		return false;
}

//**********************************************************************************************************
// compares two strings from the text
// Params:
// new_seq_start_pos - starting position of the first string
// len - number of characters to compare
// old_seq_startpos - starting position of the second string
//**********************************************************************************************************
bool CMappingCore::compare_str_to_str_in_text(uint32_t new_seq_start_pos, uint32_t len, uint32_t old_seq_start_pos)
{
	uint32_t i = 0;
	uint32_t len_in_full_bytes = len / 2;
	uchar_t *ptr_1, *ptr_2;

	if ((new_seq_start_pos ^ old_seq_start_pos) & 1)
	{
		//--------------------------------------------------------
		//new_seq_start_pos and old_seq_start_pos differ in parity
		//--------------------------------------------------------

		if (new_seq_start_pos & 1)
		{
			//new_seq_start_pos is odd, ptr_1 points it
			ptr_1 = ref_ptr + new_seq_start_pos / 2;
			ptr_2 = ref_ptr + old_seq_start_pos / 2;
		}
		else
		{
			//old_seq_start_pos is odd, ptr_1 point it
			ptr_1 = ref_ptr + old_seq_start_pos / 2;
			ptr_2 = ref_ptr + new_seq_start_pos / 2;
		}

		uint32_t buf = ptr_1[i];
		for (; i < len_in_full_bytes; i++)
		{
			buf = (buf << 8) + ptr_1[i + 1];
			if (((buf >> 4) & 0xff) != ptr_2[i])
				return false;
		}
		if (len & 1)
			//len is odd, so one more pair is to be checked
			return LO_NIBBLE(ptr_1[len_in_full_bytes]) == HI_NIBBLE(ptr_2[len_in_full_bytes]);
	}
	else
	{
		//--------------------------------------------------------
		//new_seq_start_pos and old_seq_start_pos are if the same parity
		//--------------------------------------------------------
		ptr_1 = ref_ptr + new_seq_start_pos / 2;
		ptr_2 = ref_ptr + old_seq_start_pos / 2;

		switch (new_seq_start_pos & 1)
		{
		case true:
			//--------------------------------------------------------
			//new_seq_start_pos and old_seq_start_pos are both odd
			//--------------------------------------------------------

			if (LO_NIBBLE(ptr_1[i]) != LO_NIBBLE(ptr_2[i])) //there is "a tail" at the front
				return false;
			i++;
			if (MEMCMP(ptr_1 + i, ptr_2 + i, len_in_full_bytes - i))
				return false;

			if (len & 1)	//len is odd
				return ptr_1[len_in_full_bytes] == ptr_2[len_in_full_bytes];
			else
				return HI_NIBBLE(ptr_1[len_in_full_bytes]) == HI_NIBBLE(ptr_2[len_in_full_bytes]);
			break;

		case false:
			//--------------------------------------------------------
			//new_seq_start_pos and old_seq_start_pos are both even
			//--------------------------------------------------------

			if (MEMCMP(ptr_1, ptr_2, len_in_full_bytes))
				return false;

			if (len & 1)	//len is odd
				return HI_NIBBLE(ptr_1[len_in_full_bytes]) == HI_NIBBLE(ptr_2[len_in_full_bytes]);
			break;
		} //end_of_switch
	} //end_of_else

	return true;
}

//**********************************************************************************************************
//checking identity at length "min_read_length" 
//"pattern_len" parameter is needed for correct calculations for rc
//**********************************************************************************************************
int32_t CMappingCore::check_basic_part_of_diff_len_reads(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
	uint32_t position_in_text, bool dir_gen_flag)
{
	int32_t cmp;

	//-------------------------------DIR GENOME ----------------------------------------
	if (dir_gen_flag)
	{
		switch (position_in_text & 1)
		{
		case 1:				//we read from an odd position in the text
			cmp = compare_dir_str_with_gen(ref_ptr + position_in_text / 2, pattern_sft_ptr, min_read_len, true);
			break;

		case 0:				//we read from an even position in the text
			cmp = compare_dir_str_with_gen(ref_ptr + position_in_text / 2, pattern_ptr, min_read_len, false);
			break;
		} //end_of_switch
	}
	else
		//--------------------------------RC GENOME ----------------------------------------
	{
		uint32_t right_end_index_in_dir_gen;				//starting point (when reading backward) in the reference 	
		right_end_index_in_dir_gen = ref_size - position_in_text - 1;

		uint32_t matching_part_right_point = pattern_len - 1;	//starting point (when reading backward) in the pattern

		switch (right_end_index_in_dir_gen & 1)
		{
		case 0:				                    //we read from an even position in the text (when reading backward)
			if (matching_part_right_point & 1)	//The odd number of the character being matched in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1) + 1, min_read_len, false);

			else                                 //Even number of the character being matched in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), min_read_len, false);
			break;

		case 1:				                    //we read from an odd position in the text (when reading backward)
			if (matching_part_right_point & 1)	//The odd number of the character being matched in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), min_read_len, true);

			else                                 //Even number of the character being matched in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1), min_read_len, true);
			break;
		} //end of switch;
	}
	//--------------------------------------------------------------------------------

	return cmp;
}

//**********************************************************************************************************
//checking identity at length "predef_match_length", placed after CUR_OFFSET 
//"pattern_len" parameter is needed for correct calculations for rc
//**********************************************************************************************************
int32_t CMappingCore::check_matching_part(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
	uint32_t position_in_text, bool dir_gen_flag)
{
	int32_t cmp;

	//-------------------------------DIR GENOME ----------------------------------------
	if (dir_gen_flag)
	{
		switch (position_in_text & 1)
		{
		case 1:					//starting with an odd position in the text
			if (CUR_OFFSET & 1)	//starting with an odd position in the pattern
				cmp = compare_dir_str_with_gen(ref_ptr + position_in_text / 2, pattern_ptr + CUR_OFFSET / 2, predef_match_length, true);

			else				//starting with an even position in the pattern
				cmp = compare_dir_str_with_gen(ref_ptr + position_in_text / 2, pattern_sft_ptr + CUR_OFFSET / 2, predef_match_length, true);
			break;

		case 0:				    //starting with an even position in the text
			if (CUR_OFFSET & 1)	//starting with an odd position in the pattern
				cmp = compare_dir_str_with_gen(ref_ptr + position_in_text / 2, pattern_sft_ptr + CUR_OFFSET / 2 + 1, predef_match_length, false);
			else				//starting with an even position in the pattern
				cmp = compare_dir_str_with_gen(ref_ptr + position_in_text / 2, pattern_ptr + CUR_OFFSET / 2, predef_match_length, false);
			break;
		} //end_of_switch
	}
	else
		//--------------------------------RC GENOME ----------------------------------------
	{
		uint32_t right_end_index_in_dir_gen;	//starting point (when reading backward) in the reference 	
		right_end_index_in_dir_gen = ref_size - position_in_text - 1;

		uint32_t matching_part_right_point = pattern_len - CUR_OFFSET - 1;	//starting point (when reading backward) in the pattern

		switch (right_end_index_in_dir_gen & 1)
		{
		case 0:				                    //starting with an even position in the text (when reading backward)
			if (matching_part_right_point & 1)	//starting with an odd position in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1) + 1, predef_match_length, false);

			else                                //starting with an even position in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), predef_match_length, false);
			break;

		case 1:				                    //starting with an odd position in the text (when reading backward)
			if (matching_part_right_point & 1)	//starting with an odd position in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_ptr + (matching_part_right_point >> 1), predef_match_length, true);

			else                                //starting with an even position in the pattern
				cmp = compare_rc_str_backward(ref_ptr + (right_end_index_in_dir_gen >> 1), pattern_sft_ptr + (matching_part_right_point >> 1), predef_match_length, true);
			break;
		} //end of switch;
	}
	//--------------------------------------------------------------------------------

	return cmp;
}

// EOF 