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
	prev_dp_ptr = new uint32_t[max_query_len + 2*max_ed + 3];
	curr_dp_ptr = new uint32_t[max_query_len + 2*max_ed + 3];
	
	genome_prefetch = new uchar_t[max(max_query_len, max_text_len) + 2*max_ed + 5];

	// Structures for Myers bit-par algorithm - 64-bit variant
#ifdef LEV_MYERS_ASSERTIONS
	bp_n_words = (max_query_len + 64) / 64;

	raw_bp_PM = new uint64_t[bp_n_words * 16];		// 15 - largest alphabet code
	bp_PM = new uint64_t*[16];

	for(uint32_t i = 0; i < 16; ++i)
		bp_PM[i] = &raw_bp_PM[i*bp_n_words];

	bp_scores = new uint32_t[max_text_len+1];

	bp_M     = new bp_t*[max_text_len+1];
	bp_raw_M = new bp_t[(max_text_len+1) * bp_n_words];

	alloc_text_M = (max_text_len+1) * bp_n_words;

	for(uint32_t i = 0; i < max_text_len+1; ++i)
		bp_M[i] = &bp_raw_M[i * bp_n_words];

	for(uint32_t i = 0; i < bp_n_words; ++i)
	{
		bp_M[0][i].VP = ~(0ull);
		bp_M[0][i].VN = 0;
	}

	// Structures for Myers bit-par algorithm - 128-bit variant (SSE2)
	bp128_n_words = (max_query_len + 128) / 128;
	bp128_n_chars = (bp128_n_words * 2 + 15) / 16;

	bp128_raw_PM = new Vec2uq[bp128_n_words * 16 + 8];		// 15 - largest alphabet code
	bp128_PM = new Vec2uq*[16];

	// Align to 128B boundary
	Vec2uq *tmp = bp128_raw_PM;
//	while(((uint64_t) tmp) % 128)
//		++tmp;
	for(uint32_t i = 0; i < 16; ++i)
		bp128_PM[i] = &tmp[i*bp128_n_words];

	bp128_M     = new bp128_t*[max_text_len+1 + bp128_n_words];
	bp128_raw_M = new bp128_t[(max_text_len+1 + bp128_n_words) * bp128_n_words];

//	alloc_text_M = (max_text_len+1) * bp128_n_words;

	for(uint32_t i = 0; i < max_text_len+1; ++i)
		bp128_M[i] = &bp128_raw_M[i * bp128_n_words];

	for(uint32_t i = 0; i < bp128_n_words; ++i)
	{
		bp128_M[0][i].VP = ~(0ull);
		bp128_M[0][i].VN = 0;

		bp128_M[0][i].HN = 0;
		bp128_M[0][i].HP = 0;
		bp128_M[0][i].D0 = ~(0ull);
	}

	bp128_chars = new Vec16c[bp128_n_chars];		
#endif
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

#ifdef LEV_MYERS_ASSERTIONS
	delete[] bp_scores;

	delete[] raw_bp_PM;
	delete[] bp_PM;
	delete[] bp_M;
	delete[] bp_raw_M;

	delete[] bp128_raw_PM;
	delete[] bp128_PM;
	delete[] bp128_M;
	delete[] bp128_raw_M;
	delete[] bp128_chars;
#endif
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
//				errors_in_prefix = max_mismatches;
//				break;
			}
			// More than one diagonal to compute
			else
			{
				pattern_character = GET_CHAR_FROM_PATTERN(strlenb - 1 - i, pattern_ptr);

		/*		curr_dp_ptr[left_side] = (pattern_character == genome_prefetch[left_side]) ? prev_dp_ptr[left_side-1] : 1+min(prev_dp_ptr[left_side], prev_dp_ptr[left_side-1]);
		
				for(int32_t j = left_side + 1; j < right_side; ++j)
					curr_dp_ptr[j] = (pattern_character == genome_prefetch[j]) ? prev_dp_ptr[j-1] : 1+min(min(curr_dp_ptr[j-1], prev_dp_ptr[j]), prev_dp_ptr[j-1]);

				genome_prefetch[right_side] = GET_CHAR_FROM_GENOME(genome_offset - right_side);
				curr_dp_ptr[right_side] = (pattern_character == genome_prefetch[right_side]) ? prev_dp_ptr[right_side-1] : 1+min(curr_dp_ptr[right_side-1], prev_dp_ptr[right_side-1]);
				*/

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

//		int32_t strlena = pattern_len + max_mismatches_in_suffix - fixed_right;
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

/*				curr_dp_ptr[left_side] = (pattern_character == genome_prefetch[left_side]) ? prev_dp_ptr[left_side-1] : 1+min(prev_dp_ptr[left_side], prev_dp_ptr[left_side-1]);
		
				for(int32_t j = left_side + 1; j < right_side; ++j)
					curr_dp_ptr[j] = (pattern_character == genome_prefetch[j]) ? prev_dp_ptr[j-1] : 1+min(min(curr_dp_ptr[j-1], prev_dp_ptr[j]), prev_dp_ptr[j-1]);

				genome_prefetch[right_side] = GET_CHAR_FROM_GENOME(genome_offset + right_side);
				curr_dp_ptr[right_side] = (pattern_character == genome_prefetch[right_side]) ? prev_dp_ptr[right_side-1] : 1+min(curr_dp_ptr[right_side-1], prev_dp_ptr[right_side-1]);
*/
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

#ifdef LEV_MYERS_ASSERTIONS
// ************************************************************************************
bool CEditDist::LevMyers_PP(uchar_t *seq, uint32_t _seq_len, genome_t orientation)
{
	seq_len = _seq_len;

	// Preprocessing
	fill_n(raw_bp_PM, bp_n_words * 16, 0);

	if(orientation == genome_t::direct)
		for(uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = GET_CHAR_FROM_PATTERN(i, seq);
			bp_PM[c][i / 64] |= 1ull << (i % 64);
		}
	else
		for(uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = rev_comp_code[GET_CHAR_FROM_PATTERN(seq_len-i-1, seq)];
			bp_PM[c][i / 64] |= 1ull << (i % 64);
		}

	return true;
}

// ************************************************************************************
bool CEditDist::LevMyers_PP_rawseq(uchar_t *seq, uint32_t _seq_len, genome_t orientation)
{
	seq_len = _seq_len;

	// Preprocessing
	fill_n(raw_bp_PM, bp_n_words * 16, 0);

	if (orientation == genome_t::direct)
		for (uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = raw_code[seq[i]];
			bp_PM[c][i / 64] |= 1ull << (i % 64);
		}
	else
		for (uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = raw_rev_comp_code[seq[seq_len - i - 1]];
			bp_PM[c][i / 64] |= 1ull << (i % 64);
		}

	return true;
}

// ************************************************************************************
bool CEditDist::LevMyers128_PP(uchar_t *seq, uint32_t _seq_len, genome_t orientation)
{
	return LevMyers_PP(seq, _seq_len, orientation);

	seq_len = _seq_len;

	// Preprocessing
	fill_n(bp128_raw_PM, bp128_n_words * 16, Vec2uq(0, 0));

	if(orientation == genome_t::direct)
		for(uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = GET_CHAR_FROM_PATTERN(i, seq);
			bp128_PM[c][i / 128].set_bit(i % 128, 1);
		}
	else
		for(uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = rev_comp_code[GET_CHAR_FROM_PATTERN(seq_len-i-1, seq)];
			bp128_PM[c][i / 128].set_bit(i % 128, 1);
		}

	return true;
}

// ************************************************************************************
// Myers bit-par algorithm 
bool CEditDist::LevMyers(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance)
{
	uint32_t j, r;

	uint32_t left_end_index_in_gen = ref_pos;

//	bp_scores[0] = seq_len;

	// Searching
	uint64_t X;
	uint32_t red_m = (seq_len - 1) % 64;
	uint64_t mask_m = 1ull << red_m;

	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;

	bp_t *curr_bp_M = bp_M[1];
	bp_t *prev_bp_M = bp_M[0];

	if(max_distance_in_ref == 0)
		max_distance_in_ref = max_text_len;

	for(j = 1; j <= max_distance_in_ref; ++j)
	{
		uint32_t c = GET_CHAR_FROM_GENOME(j-1);

		if(curr_bp_M - bp_raw_M >= alloc_text_M)
		{
			cerr << "To far in bp_M!\n";
			cerr << curr_bp_M - bp_raw_M << "   " << alloc_text_M << "\n";
			exit(1);
		}

//		continue;
		// r == 0
//		X = bp_PM[text[j-1]][0];
		X = bp_PM[c][0];
		curr_bp_M->D0 = (((X & prev_bp_M->VP) + prev_bp_M->VP) ^ prev_bp_M->VP) | X | prev_bp_M->VN;
		curr_bp_M->HP = prev_bp_M->VN | ~(curr_bp_M->D0 | prev_bp_M->VP);
		curr_bp_M->HN = curr_bp_M->D0 & prev_bp_M->VP;
		X = curr_bp_M->HP << 1;
		curr_bp_M->VP = (curr_bp_M->HN << 1) | ~(curr_bp_M->D0 | X);
		curr_bp_M->VN = curr_bp_M->D0 & X;

		// r == 1
/*		if(bp_n_words > 0)
		{
			X = bp_PM[text[j-1]][1];
			if(bp_M[j][0].HN & (1ull << 63))
				X |= 1;
			bp_M[j][1].D0 = (((X & bp_M[j-1][1].VP) + bp_M[j-1][1].VP) ^ bp_M[j-1][1].VP) | X | bp_M[j-1][1].VN;
			bp_M[j][1].HP = bp_M[j-1][1].VN | ~(bp_M[j][1].D0 | bp_M[j-1][1].VP);
			bp_M[j][1].HN = bp_M[j][1].D0 & bp_M[j-1][1].VP;
			X = bp_M[j][1].HP << 1;
			if(bp_M[j][0].HP & (1ull << 63))
				X |= 1;
			bp_M[j][1].VP = (bp_M[j][1].HN << 1) | ~(bp_M[j][1].D0 | X);
			if(bp_M[j][0].HN & (1ull << 63))
				bp_M[j][1].VP |= 1;
			bp_M[j][1].VN = bp_M[j][1].D0 & X;
		}*/

		// r > 1
		for(r = 1; r < bp_n_words; ++r)
		{
			curr_bp_M++;
			prev_bp_M++;
			
			if(curr_bp_M - bp_raw_M >= alloc_text_M)
			{
				cerr << "To far in bp_M!\n";
				cerr << curr_bp_M - bp_raw_M << "   " << alloc_text_M << "\n";
				exit(1);
			}

			X = bp_PM[c][r];
			X |= (curr_bp_M-1)->HN >> 63;
			curr_bp_M->D0 = (((X & prev_bp_M->VP) + prev_bp_M->VP) ^ prev_bp_M->VP) | X | prev_bp_M->VN;
			curr_bp_M->HP = prev_bp_M->VN | ~(curr_bp_M->D0 | prev_bp_M->VP);
			curr_bp_M->HN = curr_bp_M->D0 & prev_bp_M->VP;
			X = curr_bp_M->HP << 1;
			X |= (curr_bp_M-1)->HP >> 63;
			curr_bp_M->VP = (curr_bp_M->HN << 1) | ~(curr_bp_M->D0 | X);
			curr_bp_M->VP |= (curr_bp_M-1)->HN >> 63;
			curr_bp_M->VN = curr_bp_M->D0 & X;
		}

		if(curr_bp_M->HP & mask_m)
			curr_ed++;
		else if(curr_bp_M->HN & mask_m)
		{
			if(--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j;
			}
		}

		curr_bp_M++;
		prev_bp_M++;
	}

	edit_distance = min_ed;

	if (edit_distance > max_mate_edit_distance)
	{
		pos = 0;
		return false;

	}
	pos = min_ed_pos;

	return true;
}

// ************************************************************************************
// Myers bit-par algorithm - 128-bit variant (SSE2)
bool CEditDist::LevMyers128(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance)
{
	uint32_t j;

	bp_scores[0] = seq_len;

	// Searching
	Vec2uq X;
	uint32_t red_m = (seq_len - 1) % 128;
//	uint64_t mask_m = 1ull << red_m;

	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;// - bp128_n_words + 1;

	bp128_t *curr_bp_M = bp128_M[1];
	bp128_t *prev_bp_M = bp128_M[0];

	if(max_distance_in_ref == 0)
		max_distance_in_ref = max_text_len;

	for(uint32_t i = 0; i < bp128_n_words; ++i)
		bp128_chars[i] = 7;			// 7 is a symbols for which all bit vectors in bp128_PM are 0s

	Vec2uq flows;

	// !!! Caly kod tej funkcji jest przy zalozeniu, ze seq_len <= 128
//	uint32_t chars[2] = {7, 7};

	// Genome prefetch
	uchar_t *ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;
	
	*ptr++ = 7;			// No symbol
	if(ref_pos & 1)
		*ptr++ = *genome_ptr++ & 0x0f;
	for(uint32_t i = 0; i < text_len_div2; ++i)
	{
		*ptr++ = *genome_ptr >> 4;
		*ptr++ = *genome_ptr++ & 0x0f;
	}

	ptr = genome_prefetch;

	// j == 1
	X = Vec2uq(bp_PM[*(ptr+1)][0], bp_PM[*ptr][1]);
	++ptr; 

	Vec2uq tmp = permute2uq<-1, 0>((*prev_bp_M).HN);
	X |= tmp >> 63;

	curr_bp_M->D0 = (((X & prev_bp_M->VP) + prev_bp_M->VP) ^ prev_bp_M->VP) | X | prev_bp_M->VN;
	curr_bp_M->HP = prev_bp_M->VN | ~(curr_bp_M->D0 | prev_bp_M->VP);
	curr_bp_M->HN = curr_bp_M->D0 & prev_bp_M->VP;
	X = curr_bp_M->HP << 1;

	tmp = permute2uq<-1, 0>((*prev_bp_M).HP[0]);
	X |= tmp >> 63;

	curr_bp_M->VP = (curr_bp_M->HN << 1) | ~(curr_bp_M->D0 | X);
	tmp = permute2uq<-1, 0>((*prev_bp_M).HN[0]);
	curr_bp_M->VP |= tmp >> 63;

	curr_bp_M->VN = curr_bp_M->D0 & X;

	bp128_M[1][0].VP.insert(1, ~(0ull));
	bp128_M[1][0].VN.insert(1, 0);

	curr_bp_M++;
	prev_bp_M++;

	// j > 1
	for(j = 2; j <= max_distance_in_ref + 1; ++j)
	{
		// r == 0
		X = Vec2uq(bp_PM[*(ptr+1)][0], bp_PM[*ptr][1]);
		++ptr; 

		Vec2uq tmp = permute2uq<-1, 0>((*prev_bp_M).HN);
		X |= tmp >> 63;

		curr_bp_M->D0 = (((X & prev_bp_M->VP) + prev_bp_M->VP) ^ prev_bp_M->VP) | X | prev_bp_M->VN;
		curr_bp_M->HP = prev_bp_M->VN | ~(curr_bp_M->D0 | prev_bp_M->VP);
		curr_bp_M->HN = curr_bp_M->D0 & prev_bp_M->VP;
		X = curr_bp_M->HP << 1;

		tmp = permute2uq<-1, 0>((*prev_bp_M).HP[0]);
		X |= tmp >> 63;

		curr_bp_M->VP = (curr_bp_M->HN << 1) | ~(curr_bp_M->D0 | X);
		tmp = permute2uq<-1, 0>((*prev_bp_M).HN[0]);
		curr_bp_M->VP |= tmp >> 63;

		curr_bp_M->VN = curr_bp_M->D0 & X;

		if(curr_bp_M->HN.get_bit(red_m))
		{
			if(--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j;
			}
		}
		else
		{
			if(curr_bp_M->HP.get_bit(red_m))
				curr_ed++;
			if(curr_ed > max_distance_in_ref + 1 - j + max_mate_edit_distance)
				break;				
		}

		curr_bp_M++;
		prev_bp_M++;
	}

	edit_distance = min_ed;

	if(edit_distance > max_mate_edit_distance)
		return false;

	pos = min_ed_pos;

	return true;
}
#endif

// EOF
