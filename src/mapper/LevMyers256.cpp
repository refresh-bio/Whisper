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


#include "LevMyers256.h"
#include "vector_utils.h"
#include "../common/types.h"

#include <iostream>
#include <algorithm>
#include <fstream>

#undef min
#undef max

using namespace std;

#define HI_NIBBLE(x)		((x) >> 4)
#define LO_NIBBLE(x)		((x) & 0x0f)

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) ((i) & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

// ************************************************************************************
LevMyers256::LevMyers256(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed)
	: LevMyers64(_max_query_len, _max_text_len, _max_ed, 4),
	genome_prefetch(nullptr),
	bp256_raw_ptr_M(nullptr),
	bp256_M(nullptr)
{
	reallocBuffers(_max_query_len, _max_text_len, 4);
}

// ************************************************************************************
bool LevMyers256::dynamicProgramming(
	ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance) 
{
	uint32_t j;

	// Searching
	uint32_t red_m = (seq_len - 1) % 256;
	Vec4uq red_m_mask = 0;
	red_m_mask.set_bit(red_m, 1);

	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;// - bp128_n_words + 1;

	if (max_distance_in_ref == 0) {
		max_distance_in_ref = max_text_len;
		throw std::runtime_error("max_distance_in_ref == 0");
	}

	if (max_distance_in_ref > max_text_len) {
		reallocBuffers(max_query_len, max_distance_in_ref, 4);
		throw std::runtime_error("reallocation needed");
	}

	bp256_t *curr_bp_M = bp256_M[1];

	// Assumption: seq_len <= 256
									// Genome prefetch
	uchar_t *ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;

	*ptr++ = 7;			// No symbol
	*ptr++ = 7;			// No symbol
	*ptr++ = 7;			// No symbol
	if (ref_pos & 1)
		*ptr++ = *genome_ptr++ & 0x0f;
	for (uint32_t i = 0; i < text_len_div2; ++i)
	{
		*ptr++ = *genome_ptr >> 4;
		*ptr++ = *genome_ptr++ & 0x0f;
	}

	ptr = genome_prefetch;

	simd256_t prev_D0 = bp256_M[0][0].D0;
	simd256_t prev_HN = bp256_M[0][0].HN;
	simd256_t prev_HP = bp256_M[0][0].HP;
	simd256_t prev_VN = bp256_M[0][0].VN;
	simd256_t prev_VP = bp256_M[0][0].VP;

	simd256_t curr_D0;
	simd256_t curr_HN;
	simd256_t curr_HP;
	simd256_t curr_VP;
	simd256_t curr_VN;

 	for (j = 1; j <= max_distance_in_ref + 1; ++j)
	{
		// r == 0
		simd256_t X(bp_PM[*(ptr + 3)][0], bp_PM[*(ptr + 2)][1], bp_PM[*(ptr + 1)][2], bp_PM[*ptr][3]);
		++ptr;

		simd256_t tmp = permute4uq<-1, 0, 1, 2>(prev_HN);
		X |= tmp >> 63;

		curr_D0 = (((X & prev_VP) + prev_VP) ^ prev_VP) | X | prev_VN;
		curr_HP = prev_VN | ~(curr_D0 | prev_VP);
		curr_HN = curr_D0 & prev_VP;
		X = curr_HP << 1;

		tmp = permute4uq<-1, 0, 1, 2>(prev_HP);
		X |= tmp >> 63;
	
		curr_VP = (curr_HN << 1) | ~(curr_D0 | X);
		tmp = permute4uq<-1, 0, 1, 2>(prev_HN);
		curr_VP |= tmp >> 63;

		curr_VN = curr_D0 & X;

		// store values in array
		curr_bp_M->D0 = curr_D0;
		curr_bp_M->HN = curr_HN;
		curr_bp_M->HP = curr_HP;
		curr_bp_M->VN = curr_VN;
		curr_bp_M->VP = curr_VP;

		prev_D0 = curr_D0;
		prev_HN = curr_HN;
		prev_HP = curr_HP;
		prev_VN = curr_VN;
		prev_VP = curr_VP;

		if (curr_HN.get_bit(red_m))
	//	if (horizontal_or(curr_bp_M->HN & red_m_mask != 0))
		{
			if (--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j - (seq_len / 64); 
			}
			else if (curr_ed == min_ed)
				if (min_ed_pos + 1 == j - (seq_len / 64))
					min_ed_pos = j - (seq_len / 64); 
		}
		else
		{
			if (curr_HP.get_bit(red_m))
		//	if (horizontal_or(curr_bp_M->HP & red_m_mask != 0))
				curr_ed++;
			else if (curr_ed == min_ed)
				if(min_ed_pos + 1 == j - (seq_len / 64))
					min_ed_pos = j - (seq_len / 64); 

			if (curr_ed > max_distance_in_ref + 1 - j + max_mate_edit_distance)
				break;
		}

		curr_bp_M++;
	}

	edit_distance = min_ed;

	if (edit_distance > max_mate_edit_distance) {
		pos = 0;
		return false;
	}

	pos = min_ed_pos;

	return true;
}

// ************************************************************************************
ref_pos_t LevMyers256::getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
	ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) 
{
	uint32_t read_pos = seq_len;
	uint32_t cigar_pos = 0;
	
	const char *decode = "ACGTNNXX";

	uint32_t h_val, v_val, d_val, c_val;

	c_val = edit_dist;
	affine_score = 0;
	bool inIndel = false;
	num_events = 0;
	ext_cigar[0] = 0;
	string_reader_t quality_reader(quality, seq_len, dir);

	while (read_pos && pos) {
		if (tmp_read_sequence[read_pos] == tmp_ref_sequence[pos - 1]) {
			// match
			inIndel = false;
			affine_score += scoring.match; //quality_reader[read_pos];
			--read_pos;
			--pos;
			ext_cigar[cigar_pos++] = '.';
		}
		else {
			uint32_t word_no = (read_pos - 1) / 64;
			uint64_t word_mask = 1ull << ((read_pos - 1) % 64);

			const auto& word256 = *bp256_M[pos + word_no];

			d_val = h_val = v_val = c_val;

			if (!(word256.D0.extract(word_no) & word_mask)) { --d_val; }
			if (word256.HP.extract(word_no) & word_mask) { --h_val; }
			if (word256.HN.extract(word_no) & word_mask) { ++h_val; }
			if (word256.VP.extract(word_no) & word_mask) { --v_val; }
			if (word256.VN.extract(word_no) & word_mask) { ++v_val; }

			if (d_val <= h_val && d_val <= v_val) {
				// mismatch
				affine_score += scoring.mismatch; //quality_reader[read_pos];
				inIndel = false;
				++num_events;

				ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos-1]];
				--read_pos;
				--pos;
				c_val = d_val;
			}
			else if (h_val <= v_val) {
				// deletion
				if (inIndel) {
					affine_score += scoring.gap_extend;
				} else {
					affine_score += scoring.gap_open;
					++num_events;
					inIndel = true;
				}

				ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos-1]];
				ext_cigar[cigar_pos++] = '^';
				--pos;
				c_val = h_val;
			}
			else {
				// insertion
				if (inIndel) {
					affine_score += scoring.gap_extend;
				} else {
					affine_score += scoring.gap_open;
					++num_events;
					inIndel = true;
				}

				ext_cigar[cigar_pos++] = decode[tmp_read_sequence[read_pos]];
				ext_cigar[cigar_pos++] = '#';
				--read_pos;
				c_val = v_val;
			}
		}
	}

	while (read_pos)
	{
		// insertion
		if (inIndel) 
		{
			affine_score += scoring.gap_extend;
		}
		else 
		{
			affine_score += scoring.gap_open;
			++num_events;
			inIndel = true;
		}
		
		ext_cigar[cigar_pos++] = decode[tmp_read_sequence[read_pos]];
		ext_cigar[cigar_pos++] = '#';
		--read_pos;
	}

	ext_cigar[cigar_pos] = '\0';

	reverse(ext_cigar, ext_cigar + cigar_pos);
	return pos;
}

// ************************************************************************************
#ifdef _DEVELOPMENT_MODE
void LevMyers256::save(const std::string& filename) {

	std::ofstream file(filename);

	for (uint32_t i = 0; i < max_text_len; ++i) {
		file << std::dec << i << "," << code2symbol[genome_prefetch[i]];
		file << endl << "D0: " << printVec4uq(bp256_M[i][0].D0)
		 << endl << "HN: " << printVec4uq(bp256_M[i][0].HN)
		 << endl << "HP: " << printVec4uq(bp256_M[i][0].HP)
		 << endl << "VN: " << printVec4uq(bp256_M[i][0].VN)
		 << endl << "VP: " << printVec4uq(bp256_M[i][0].VP)
	/*	<< endl
		<< endl << "X0: " << printVec4uq(bp256_M[i][0].X0)
		<< endl << "X1: " << printVec4uq(bp256_M[i][0].X1)
		<< endl << "X2: " << printVec4uq(bp256_M[i][0].X2)
		<< endl << "X3: " << printVec4uq(bp256_M[i][0].X3)*/
		 << endl << endl;
	}
}
#endif

// EOF
