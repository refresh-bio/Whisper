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
template <instruction_set_t instruction_set>
LevMyers256<instruction_set>::LevMyers256(uint32_t _max_text_len, uint32_t _max_ed)
	: LevMyers64(255, _max_text_len, _max_ed, 4),
	genome_prefetch(nullptr),
	bp256_raw_ptr_M(nullptr),
	bp256_M(nullptr)
{
	reallocBuffers(max_query_len, _max_text_len, 4);
}

// ************************************************************************************
template <instruction_set_t instruction_set>
void LevMyers256<instruction_set>::reallocBuffers(uint32_t _max_query_len, uint32_t _max_text_len, int rounding)
{
	this->max_query_len = _max_query_len;
	this->max_text_len = _max_text_len,

	bp256_n_words = (max_query_len + 256) / 256;
	
	if (auto t = (bp256_t**)realloc(bp256_M, sizeof(bp256_t*) * (max_text_len + 4 + bp256_n_words))) // + 4 because 256 bit
		bp256_M = t;
	else
		exit(1);

	if (auto t = (bp256_t*)alloc_aligned(bp256_raw_ptr_M, (max_text_len + 4 + bp256_n_words) * bp256_n_words * sizeof(bp256_t), sizeof(simd256_t)))
		bp256_raw_M = t;
	else
		exit(1);

	for (uint32_t i = 0; i < max_text_len + 4; ++i)
		bp256_M[i] = &bp256_raw_M[i * bp256_n_words];

	for (uint32_t i = 0; i < bp256_n_words; ++i)
	{
		bp256_M[0][i].VP = ~(0ull);
		bp256_M[0][i].VN = 0;

		bp256_M[0][i].HN = 0;
		bp256_M[0][i].HP = 0;
		bp256_M[0][i].D0 = ~(0ull);
	}

	if (auto t = (uchar_t*)realloc(genome_prefetch, sizeof(uchar_t) * (std::max(max_query_len, max_text_len) + 8)))
		genome_prefetch = t;
	else
		exit(1);
}

// ************************************************************************************
template <instruction_set_t instruction_set>
bool LevMyers256<instruction_set>::dynamicProgramming(
	ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance) 
{
	uint32_t j;

	// Searching
	uint32_t seq_len_m1 = seq_len - 1;
	uint32_t red_m = (seq_len_m1) % 256;
	uint32_t word_num = seq_len_m1 / 64;
	uint64_t red_m_mask = 1ull << (red_m % 64);
	uint64_t red_m_word = red_m / 64;

	//int middle_point_in_ref = (max_distance_in_ref + seq_len) / 2;
	//int best_dist_from_middle = max_distance_in_ref + seq_len;		// large value

	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;// - bp128_n_words + 1;

	if (max_distance_in_ref == 0) {
		max_distance_in_ref = max_text_len;
		throw std::runtime_error("max_distance_in_ref == 0");
	}

	if (max_distance_in_ref > max_text_len) {
		throw std::runtime_error("reallocation needed: max_query_len: " + to_string(max_text_len) + " max_dist: " + to_string(max_distance_in_ref));
		reallocBuffers(max_query_len, max_distance_in_ref, 4);
	}

	bp256_t *curr_bp_M = bp256_M[1];

	// Assumption: seq_len <= 256
									// Genome prefetch
	uchar_t *ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	//uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;

	*ptr++ = 7;			// No symbol
	*ptr++ = 7;			// No symbol
	*ptr++ = 7;			// No symbol
	PrefetchDecompressed256<instruction_set>(ptr, genome_ptr, max_distance_in_ref + 3, ref_pos & 1);
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

		simd256_t tmp = permute4<-1, 0, 1, 2>(prev_HN);
		X |= tmp >> 63;

		curr_D0 = (((X & prev_VP) + prev_VP) ^ prev_VP) | X | prev_VN;
		curr_HP = prev_VN | ~(curr_D0 | prev_VP);
		curr_HN = curr_D0 & prev_VP;
		X = curr_HP << 1;

		tmp = permute4<-1, 0, 1, 2>(prev_HP);
		X |= tmp >> 63;
	
		curr_VP = (curr_HN << 1) | ~(curr_D0 | X);
		tmp = permute4<-1, 0, 1, 2>(prev_HN);
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

//		if (curr_HN.get_bit(red_m))
		if (curr_HN[red_m_word] & red_m_mask)
		{
			if (--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j - word_num;
//				best_dist_from_middle = abs((int) j - middle_point_in_ref);
			}
			else if (curr_ed == min_ed)
			{
/*				if (min_ed_pos + 1 == j - word_num)
				{
					best_dist_from_middle = dist_from_middle;
					min_ed_pos = j - word_num;
				}*/

				bool last_symbol_match = *(ptr + (3 - (int)word_num) - 1) == last_symbol;

				auto D0_bit = curr_D0[red_m_word] & red_m_mask;
//				if((curr_D0.get_bit(red_m) && last_symbol_match) || !curr_D0.get_bit(red_m))
				if((D0_bit && last_symbol_match) || !D0_bit)
					min_ed_pos = j - word_num;

/*				int dist_from_middle = abs((int)j - middle_point_in_ref);
				if (min_ed_pos + 1 == j - word_num && dist_from_middle < best_dist_from_middle)
				{
					best_dist_from_middle = dist_from_middle;
					min_ed_pos = j - word_num;
				}*/
			}
		}
		else
		{
//			if (curr_HP.get_bit(red_m))
			if (curr_HP[red_m_word] & red_m_mask)
				curr_ed++;
			else if (curr_ed == min_ed)
			{
				bool last_symbol_match = *(ptr + (3 - (int) word_num) - 1) == last_symbol;

//				if ((curr_D0.get_bit(red_m) && last_symbol_match) || !curr_D0.get_bit(red_m))
				auto D0_bit = curr_D0[red_m_word] & red_m_mask;
				if ((D0_bit && last_symbol_match) || !D0_bit)
					min_ed_pos = j - word_num;


/*				int dist_from_middle = abs((int)j - middle_point_in_ref);
				if (min_ed_pos + 1 == j - word_num && dist_from_middle < best_dist_from_middle)
				{
					best_dist_from_middle = dist_from_middle;
					min_ed_pos = j - word_num;
				}*/

//				if (min_ed_pos + 1 == j - word_num)
//					min_ed_pos = j - word_num;
			}

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
template <instruction_set_t instruction_set>
ref_pos_t LevMyers256<instruction_set>::getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
	ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) 
{
	uint32_t read_pos = seq_len;
	uint32_t cigar_pos = 0;
	
	const char *decode = "ACGTNNXX";

	uint32_t h_val, v_val, d_val, c_val;
	matrix_dir_t prev_dir = matrix_dir_t::D;
	matrix_dir_t curr_dir;

	c_val = edit_dist;
	affine_score = 0;
	bool inIndel = false;
	num_events = 0;
	ext_cigar[0] = 0;
	string_reader_t quality_reader(quality, seq_len, dir);

	while (read_pos && pos) {
		uint32_t word_no = (read_pos - 1) / 64;
		uint64_t word_mask = 1ull << ((read_pos - 1) % 64);

		const auto& word256 = *bp256_M[pos + word_no];

		d_val = h_val = v_val = c_val;

		if (!(word256.D0.extract(word_no) & word_mask)) { --d_val; }
		if (word256.HP.extract(word_no) & word_mask) { --h_val; }
		if (word256.HN.extract(word_no) & word_mask) { ++h_val; }
		if (word256.VP.extract(word_no) & word_mask) { --v_val; }
		if (word256.VN.extract(word_no) & word_mask) { ++v_val; }

		int match_cost = (tmp_read_sequence[read_pos] != tmp_ref_sequence[pos - 1]);

		bool can_come_from_H = (h_val + 1 == c_val);
		bool can_come_from_V = (v_val + 1 == c_val);
		bool can_come_from_D = (d_val + match_cost == c_val);

		if (can_come_from_D && (
			(prev_dir == matrix_dir_t::D) ||
			(prev_dir == matrix_dir_t::H && !can_come_from_H) ||
			(prev_dir == matrix_dir_t::V && !can_come_from_V)))
			curr_dir = matrix_dir_t::D;
		else if (can_come_from_H && (
			(prev_dir == matrix_dir_t::H) ||
			(prev_dir == matrix_dir_t::D && !can_come_from_D) ||
			(prev_dir == matrix_dir_t::V && !can_come_from_V)))
			curr_dir = matrix_dir_t::H;
		else
			curr_dir = matrix_dir_t::V;

		if (curr_dir == matrix_dir_t::D)
		{
			if (tmp_read_sequence[read_pos] == tmp_ref_sequence[pos - 1]) {
				// match
				affine_score += scoring.match; //quality_reader[read_pos];
				ext_cigar[cigar_pos++] = '.';
			}
			else {
				// mismatch
				affine_score += scoring.mismatch; //quality_reader[read_pos];
				ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos - 1]];
			}

			inIndel = false;
			--read_pos;
			--pos;
			c_val = d_val;
		}
		else if (curr_dir == matrix_dir_t::H)
		{
			// deletion
			if (inIndel) {
				affine_score += scoring.gap_del_extend;
			}
			else {
				affine_score += scoring.gap_del_open;
				++num_events;
				inIndel = true;
			}

			ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos - 1]];
			ext_cigar[cigar_pos++] = '^';
			--pos;
			c_val = h_val;
		}
		else
		{
			// insertion
			if (inIndel) {
				affine_score += scoring.gap_ins_extend;
			}
			else {
				affine_score += scoring.gap_ins_open;
				++num_events;
				inIndel = true;
			}

			ext_cigar[cigar_pos++] = decode[tmp_read_sequence[read_pos]];
			ext_cigar[cigar_pos++] = '#';
			--read_pos;
			c_val = v_val;
		}

		prev_dir = curr_dir;
	}

	while (read_pos)
	{
		// insertion
		if (inIndel) 
		{
			affine_score += scoring.gap_ins_extend;
		}
		else 
		{
			affine_score += scoring.gap_ins_open;
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
template <instruction_set_t instruction_set>
void LevMyers256<instruction_set>::save(const std::string& filename) {


	std::ofstream file(filename);

	for (uint32_t i = 0; i < max_text_len; ++i) {
		file << std::dec << i << "," << code2symbol[genome_prefetch[i]];
		file << endl << "D0: " << printVec4uq(bp256_M[i][0].D0)
		 << endl << "HN: " << printVec4uq(bp256_M[i][0].HN)
		 << endl << "HP: " << printVec4uq(bp256_M[i][0].HP)
		 << endl << "VN: " << printVec4uq(bp256_M[i][0].VN)
		 << endl << "VP: " << printVec4uq(bp256_M[i][0].VP)
		 << endl << endl;
	}
}
#endif

// EOF
