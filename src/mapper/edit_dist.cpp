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


#include "edit_dist.h"
#include "../common/utils.h"

#include <iostream>
#include <algorithm>

using namespace std;

#define HI_NIBBLE(x)		((x) >> 4)
#define LO_NIBBLE(x)		((x) & 0x0f)

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) ((i) & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

// ************************************************************************************
CEditDist::CEditDist(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed)
{
	max_query_len = _max_query_len;
	max_text_len  = _max_text_len;
	max_ed        = _max_ed;

	prev_dp_ptr = nullptr;
	curr_dp_ptr = nullptr;
	genome_prefetch = nullptr;

	ref_ptr = nullptr;

	allocate();

	rev_comp_code[(int) sym_code_A] = sym_code_T;
	rev_comp_code[(int) sym_code_C] = sym_code_G;
	rev_comp_code[(int) sym_code_G] = sym_code_C;
	rev_comp_code[(int) sym_code_T] = sym_code_A;
	rev_comp_code[4]		  	    = 4;
	rev_comp_code[5]			    = 5;
	rev_comp_code[6]				= 6;
	rev_comp_code[7]				= 7;

	for (int i = 0; i < 128; ++i)
	{
		raw_code[i] = 6;
		raw_rev_comp_code[i] = 6;
	}

	raw_code['A'] = sym_code_A;
	raw_code['C'] = sym_code_C;
	raw_code['G'] = sym_code_G;
	raw_code['T'] = sym_code_T;
	raw_code['N'] = sym_code_N_read;

	raw_rev_comp_code['A'] = sym_code_T;
	raw_rev_comp_code['C'] = sym_code_G;
	raw_rev_comp_code['G'] = sym_code_C;
	raw_rev_comp_code['T'] = sym_code_A;
	raw_rev_comp_code['N'] = sym_code_N_read;
}

// ************************************************************************************
CEditDist::~CEditDist()
{
	release();
}

// ************************************************************************************
void CEditDist::allocate()
{
	// Structures for LevDiag
	prev_dp_ptr = new uint32_t[max_query_len + 2ull * max_ed + 3];
	curr_dp_ptr = new uint32_t[max_query_len + 2ull * max_ed + 3];
	
	genome_prefetch = new uchar_t[max(max_query_len, max_text_len) + 2ull * max_ed + 5];
}

// ************************************************************************************
void CEditDist::release()
{
	if(prev_dp_ptr)
		delete[] prev_dp_ptr;
	if(curr_dp_ptr)
		delete[] curr_dp_ptr;
	if(genome_prefetch)
		delete[] genome_prefetch;
}

// ************************************************************************************
void CEditDist::SetReference(uchar_t *_ref_ptr, uint32_t _ref_size, uint32_t _cur_offset)
{
	ref_ptr = _ref_ptr;
	ref_size = _ref_size;
	cur_offset = _cur_offset;
}

// ************************************************************************************
uint32_t CEditDist::LevDiag(uchar_t* pattern_ptr, uint32_t pattern_len, uint32_t position_in_text,
								uint32_t fixed_left, uint32_t fixed_right, uint32_t max_mismatches, bool dir_gen_flag)
{
	uint32_t left_end_index_in_gen; 
	uint8_t pattern_character;

	if(dir_gen_flag)
		left_end_index_in_gen = position_in_text - cur_offset - max_mismatches;
	else
		left_end_index_in_gen = ref_size - position_in_text - pattern_len + cur_offset - max_mismatches;
  
	int errors_in_prefix = 0;
	int errors_in_suffix = 0;
	
	// Testing of fixed segment
	for(uint32_t i = fixed_left; i < fixed_right; ++i)
		if (GET_CHAR_FROM_GENOME(i + max_mismatches) != GET_CHAR_FROM_PATTERN(i, pattern_ptr))
		{
			return max_mismatches + 1;

			string s_err;
			for (uint32_t j = fixed_left; j < fixed_right; ++j)
			{
				s_err.push_back("ACGTNxyz"[GET_CHAR_FROM_GENOME(j + max_mismatches)]);
				s_err.push_back("ACGTNxyz"[GET_CHAR_FROM_PATTERN(j, pattern_ptr)]);
				s_err.push_back(' ');
			}

			s_err.push_back('\n');
			if (dir_gen_flag)
				s_err.push_back('d');
			else
				s_err.push_back('r');
			for (uint32_t j = 0; j < pattern_len; ++j)
				s_err.push_back("ACGTNxyz"[GET_CHAR_FROM_PATTERN(j, pattern_ptr)]);

//			cerr << "Mismatch in fixed segment: " + s_err + "\n";

//			break;
		}
		

	// *** Processing pattern prefix 
	// Try to extend fixed region to the left
	while(fixed_left)
		if(GET_CHAR_FROM_GENOME(fixed_left + max_mismatches - 1) == GET_CHAR_FROM_PATTERN(fixed_left-1, pattern_ptr))
			--fixed_left;
		else
			break;

	if(fixed_left > 0)
	{
		int32_t strlena = fixed_left + max_mismatches;
		int32_t strlenb = fixed_left;
		uint32_t genome_offset = strlena + max_mismatches;

		int32_t left_side  = 1;
		int32_t right_side = 2 * max_mismatches + 1;

		for(uint32_t j = 0; j < max_mismatches+1; ++j)
			prev_dp_ptr[max_mismatches + j] = prev_dp_ptr[max_mismatches - j] = j;

		// Genome prefetch (only prefix)
		for(int32_t j = left_side; j < right_side; ++j)
			genome_prefetch[j] = GET_CHAR_FROM_GENOME(genome_offset - j);

		for(int32_t i = 0; i < strlenb; i++)
		{
			// Single diagonal to compute
			if(left_side == right_side)
			{
				if(prev_dp_ptr[left_side-1] > max_mismatches)
					return max_mismatches + 1;

				for(; i < strlenb; ++i)
				{
					if(GET_CHAR_FROM_PATTERN(strlenb-1-i, pattern_ptr) != GET_CHAR_FROM_GENOME(genome_offset - left_side))
						return max_mismatches + 1;
					++left_side;
				}
				--left_side;
				curr_dp_ptr[left_side] = max_mismatches;
				right_side = left_side;
			}
			// More than one diagonal to compute
			else
			{
				pattern_character = GET_CHAR_FROM_PATTERN(strlenb - 1 - i, pattern_ptr);

				uint32_t lower_left = prev_dp_ptr[left_side-1];
				for(int32_t j = left_side; j < right_side; ++j)
				{
					if(pattern_character == genome_prefetch[j])
					{
						curr_dp_ptr[j] = prev_dp_ptr[j-1];
						lower_left = min(curr_dp_ptr[j], prev_dp_ptr[j]);
					}
					else
					{
						if(prev_dp_ptr[j] > lower_left)
							curr_dp_ptr[j] = ++lower_left;
						else
						{
							lower_left = prev_dp_ptr[j];
							curr_dp_ptr[j] = lower_left+1;
						}
					}
				}

				genome_prefetch[right_side] = GET_CHAR_FROM_GENOME(genome_offset - right_side);
				if(pattern_character == genome_prefetch[right_side])
					curr_dp_ptr[right_side] = prev_dp_ptr[right_side-1];
				else
					curr_dp_ptr[right_side] = lower_left+1;

				// Narrow down the range of diagonals if possible
				if(curr_dp_ptr[left_side] > max_mismatches)
					if(++left_side > right_side)
						return max_mismatches + 1;
				if(curr_dp_ptr[right_side] > max_mismatches)
					if(left_side > --right_side)
						return max_mismatches + 1;
			}
			
			swap(curr_dp_ptr, prev_dp_ptr);
			++left_side;
			++right_side;
		}  

		if(left_side > right_side)
			return max_mismatches + 1;
			
		errors_in_prefix = max_mismatches+1;
		for(int j = left_side-1; j < right_side; ++j)
			errors_in_prefix = min(errors_in_prefix, (int32_t) prev_dp_ptr[j]);
	
		if(errors_in_prefix > (int32_t) max_mismatches)
			return max_mismatches + 1;
	}

	// *** Processing pattern suffix
	// Try to extend fixed region to the right
	while(fixed_right < pattern_len)
		if(GET_CHAR_FROM_GENOME(fixed_right + max_mismatches) == GET_CHAR_FROM_PATTERN(fixed_right, pattern_ptr))
			++fixed_right;
		else
			break;

	if(fixed_right < pattern_len)
	{
		int32_t max_mismatches_in_suffix = max_mismatches - errors_in_prefix;

		int32_t strlenb = pattern_len - fixed_right;
  
		int left_side  = 1;
		int right_side = 2 * max_mismatches_in_suffix + 1;

		for(int32_t j = 0; j < max_mismatches_in_suffix+1; ++j)
			prev_dp_ptr[max_mismatches_in_suffix + j] = prev_dp_ptr[max_mismatches_in_suffix - j] = j;

		uint32_t genome_offset = max_mismatches + fixed_right - max_mismatches_in_suffix - 1;

		// Genome prefetch (only prefix)
		for(int32_t j = left_side; j < right_side; ++j)
			genome_prefetch[j] = GET_CHAR_FROM_GENOME(genome_offset + j);

		for(int i = 0; i < strlenb; i++)
		{
			// Single diagonal to compute
			if(left_side == right_side)
			{
				if((int32_t) prev_dp_ptr[left_side-1] > max_mismatches_in_suffix)
					return max_mismatches + 1;

				for(; i < strlenb; ++i)
				{
					if(GET_CHAR_FROM_PATTERN(i + fixed_right, pattern_ptr) != GET_CHAR_FROM_GENOME(genome_offset + left_side))
						return max_mismatches + 1;
					++left_side;
				}
				--left_side;
				curr_dp_ptr[left_side] = max_mismatches_in_suffix;
				right_side = left_side;
			}
			// More than one diagonal to compute
			else
			{
				pattern_character = GET_CHAR_FROM_PATTERN(i + fixed_right, pattern_ptr);

				uint32_t lower_left = prev_dp_ptr[left_side-1];
				for(int32_t j = left_side; j < right_side; ++j)
				{
					if(pattern_character == genome_prefetch[j])
					{
						curr_dp_ptr[j] = prev_dp_ptr[j-1];
						lower_left = min(curr_dp_ptr[j], prev_dp_ptr[j]);
					}
					else
					{
						if(prev_dp_ptr[j] > lower_left)
							curr_dp_ptr[j] = ++lower_left;
						else
						{
							lower_left = prev_dp_ptr[j];
							curr_dp_ptr[j] = lower_left+1;
						}
					}
				}

				genome_prefetch[right_side] = GET_CHAR_FROM_GENOME(genome_offset + right_side);
				if(pattern_character == genome_prefetch[right_side])
					curr_dp_ptr[right_side] = prev_dp_ptr[right_side-1];
				else
					curr_dp_ptr[right_side] = lower_left+1;

				// Narrow down the range of diagonals if possible
				if((int32_t) curr_dp_ptr[left_side] > max_mismatches_in_suffix)
					if(++left_side > right_side)
						return max_mismatches + 1;
				if((int32_t) curr_dp_ptr[right_side] > max_mismatches_in_suffix)
					if(left_side > --right_side)	
						return max_mismatches + 1;
			}

			swap(curr_dp_ptr, prev_dp_ptr);
			++left_side;
			++right_side;
		}

		if(left_side > right_side)
			return max_mismatches + 1;

		errors_in_suffix = max_mismatches_in_suffix+1;
		for(int j = left_side-1; j < right_side; ++j)
			errors_in_suffix = min((int32_t) errors_in_suffix, (int32_t) prev_dp_ptr[j]);
	}

	return errors_in_prefix + errors_in_suffix;
}

// ************************************************************************************
// Version for matching with variant indel
// variant_len > 0 means insertion
// variant_len < 0 means deletion
pair<uint32_t, int32_t> CEditDist::LevDiagVariant(uchar_t* pattern_ptr, uint32_t pattern_len, uint32_t position_in_text,
	uint32_t fixed_left, uint32_t fixed_right, int variant_len, uint32_t max_mismatches, bool dir_gen_flag)
{
	uint32_t left_end_index_in_gen;
	uint8_t pattern_character;
	int32_t left_best_pos = 0;

	if (dir_gen_flag)
		left_end_index_in_gen = position_in_text - fixed_left - max_mismatches;
	else
		left_end_index_in_gen = ref_size - position_in_text - pattern_len + fixed_left - max_mismatches;

	int errors_in_prefix = 0;
	int errors_in_suffix = 0;

	// *** Processing pattern prefix 
	// Try to extend fixed region to the left
	while (fixed_left)
		if (GET_CHAR_FROM_GENOME(fixed_left + max_mismatches - 1) == GET_CHAR_FROM_PATTERN(fixed_left - 1, pattern_ptr))
			--fixed_left;
		else
			break;

	if (fixed_left > 0)
	{
		int32_t strlena = fixed_left + max_mismatches;
		int32_t strlenb = fixed_left;
		uint32_t genome_offset = strlena + max_mismatches;

		int32_t left_side = 1;
		int32_t right_side = 2 * max_mismatches + 1;

		for (uint32_t j = 0; j < max_mismatches + 1; ++j)
			prev_dp_ptr[max_mismatches + j] = prev_dp_ptr[max_mismatches - j] = j;

		// Genome prefetch (only prefix)
		for (int32_t j = left_side; j < right_side; ++j)
			genome_prefetch[j] = GET_CHAR_FROM_GENOME(genome_offset - j);

		for (int32_t i = 0; i < strlenb; i++)
		{
			// Single diagonal to compute
			if (left_side == right_side)
			{
				if (prev_dp_ptr[left_side - 1] > max_mismatches)
					return make_pair(max_mismatches + 1, 0);

				for (; i < strlenb; ++i)
				{
					if (GET_CHAR_FROM_PATTERN(strlenb - 1 - i, pattern_ptr) != GET_CHAR_FROM_GENOME(genome_offset - left_side))
						return make_pair(max_mismatches + 1, 0);
					++left_side;
				}
				--left_side;
				curr_dp_ptr[left_side] = max_mismatches;
				right_side = left_side;
				//				errors_in_prefix = max_mismatches;
				//				break;
			}
			// More than one diagonal to compute
			else
			{
				pattern_character = GET_CHAR_FROM_PATTERN(strlenb - 1 - i, pattern_ptr);

				uint32_t lower_left = prev_dp_ptr[left_side - 1];
				for (int32_t j = left_side; j < right_side; ++j)
				{
					if (pattern_character == genome_prefetch[j])
					{
						curr_dp_ptr[j] = prev_dp_ptr[j - 1];
						lower_left = min(curr_dp_ptr[j], prev_dp_ptr[j]);
					}
					else
					{
						if (prev_dp_ptr[j] > lower_left)
							curr_dp_ptr[j] = ++lower_left;
						else
						{
							lower_left = prev_dp_ptr[j];
							curr_dp_ptr[j] = lower_left + 1;
						}
					}
				}

				genome_prefetch[right_side] = GET_CHAR_FROM_GENOME(genome_offset - right_side);
				if (pattern_character == genome_prefetch[right_side])
					curr_dp_ptr[right_side] = prev_dp_ptr[right_side - 1];
				else
					curr_dp_ptr[right_side] = lower_left + 1;

				// Narrow down the range of diagonals if possible
				if (curr_dp_ptr[left_side] > max_mismatches)
					if (++left_side > right_side)
						return make_pair(max_mismatches + 1, 0);
				if (curr_dp_ptr[right_side] > max_mismatches)
					if (left_side > --right_side)
						return make_pair(max_mismatches + 1, 0);
			}

			swap(curr_dp_ptr, prev_dp_ptr);
			++left_side;
			++right_side;
		}

		if (left_side > right_side)
			return make_pair(max_mismatches + 1, 0);

		errors_in_prefix = max_mismatches + 1;
		for (int j = left_side - 1; j < right_side; ++j)
			if (errors_in_prefix > (int32_t)prev_dp_ptr[j])
			{
				errors_in_prefix = (int32_t)prev_dp_ptr[j];
				left_best_pos = j;
			}

		if (errors_in_prefix > (int32_t) max_mismatches)
			return make_pair(max_mismatches + 1, 0);
	}

	// *** Processing pattern suffix
	// Try to extend fixed region to the right
	while (fixed_right < pattern_len)
		if (GET_CHAR_FROM_GENOME(fixed_right + max_mismatches) == GET_CHAR_FROM_PATTERN(fixed_right, pattern_ptr))
			++fixed_right;
		else
			break;

	if (fixed_right < pattern_len)
	{
		int32_t max_mismatches_in_suffix = max_mismatches - errors_in_prefix;

		//		int32_t strlena = pattern_len + max_mismatches_in_suffix - fixed_right;
		int32_t strlenb = pattern_len - fixed_right;

		int left_side = 1;
		int right_side = 2 * max_mismatches_in_suffix + 1;

		for (int32_t j = 0; j < max_mismatches_in_suffix + 1; ++j)
			prev_dp_ptr[max_mismatches_in_suffix + j] = prev_dp_ptr[max_mismatches_in_suffix - j] = j;

		uint32_t genome_offset = max_mismatches + fixed_right - max_mismatches_in_suffix - 1;

		// Genome prefetch (only prefix)
		for (int32_t j = left_side; j < right_side; ++j)
			genome_prefetch[j] = GET_CHAR_FROM_GENOME(genome_offset + j);

		for (int i = 0; i < strlenb; i++)
		{
			// Single diagonal to compute
			if (left_side == right_side)
			{
				if ((int32_t)prev_dp_ptr[left_side - 1] > max_mismatches_in_suffix)
					return make_pair(max_mismatches + 1, 0);

				for (; i < strlenb; ++i)
				{
					if (GET_CHAR_FROM_PATTERN(i + fixed_right, pattern_ptr) != GET_CHAR_FROM_GENOME(genome_offset + left_side))
						return make_pair(max_mismatches + 1, 0);
					++left_side;
				}
				--left_side;
				curr_dp_ptr[left_side] = max_mismatches_in_suffix;
				right_side = left_side;
			}
			// More than one diagonal to compute
			else
			{
				pattern_character = GET_CHAR_FROM_PATTERN(i + fixed_right, pattern_ptr);

				uint32_t lower_left = prev_dp_ptr[left_side - 1];
				for (int32_t j = left_side; j < right_side; ++j)
				{
					if (pattern_character == genome_prefetch[j])
					{
						curr_dp_ptr[j] = prev_dp_ptr[j - 1];
						lower_left = min(curr_dp_ptr[j], prev_dp_ptr[j]);
					}
					else
					{
						if (prev_dp_ptr[j] > lower_left)
							curr_dp_ptr[j] = ++lower_left;
						else
						{
							lower_left = prev_dp_ptr[j];
							curr_dp_ptr[j] = lower_left + 1;
						}
					}
				}

				genome_prefetch[right_side] = GET_CHAR_FROM_GENOME(genome_offset + right_side);
				if (pattern_character == genome_prefetch[right_side])
					curr_dp_ptr[right_side] = prev_dp_ptr[right_side - 1];
				else
					curr_dp_ptr[right_side] = lower_left + 1;

				// Narrow down the range of diagonals if possible
				if ((int32_t)curr_dp_ptr[left_side] > max_mismatches_in_suffix)
					if (++left_side > right_side)
						return make_pair(max_mismatches + 1, 0);
				if ((int32_t)curr_dp_ptr[right_side] > max_mismatches_in_suffix)
					if (left_side > --right_side)
						return make_pair(max_mismatches + 1, 0);
			}

			swap(curr_dp_ptr, prev_dp_ptr);
			++left_side;
			++right_side;
		}

		if (left_side > right_side)
			return make_pair(max_mismatches + 1, 0);

		errors_in_suffix = max_mismatches_in_suffix + 1;
		for (int j = left_side - 1; j < right_side; ++j)
			errors_in_suffix = min((int32_t)errors_in_suffix, (int32_t)prev_dp_ptr[j]);
	}

	return make_pair(errors_in_prefix + errors_in_suffix, left_best_pos - (int32_t) max_mismatches);
}

// EOF
