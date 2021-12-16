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


#include "LevMyers128.h"
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

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) (i & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

// ************************************************************************************
template <instruction_set_t instruction_set>
LevMyers128<instruction_set>::LevMyers128(uint32_t _max_text_len, uint32_t _max_ed)
	:
	LevMyers64(127, _max_text_len, _max_ed, 2),
	genome_prefetch(nullptr),
	bp128_raw_ptr_M(nullptr),
	bp128_M(nullptr)
{
	reallocBuffers(max_query_len, _max_text_len, 2);
}

// ************************************************************************************
template <instruction_set_t instruction_set>
void LevMyers128<instruction_set>::reallocBuffers(uint32_t _max_query_len, uint32_t _max_text_len, int rounding)
{
	LevMyers64::reallocBuffers(_max_query_len, _max_text_len, rounding);

	bp128_n_words = (this->max_query_len + 128) / 128;

	if (auto t = (bp128_t**)realloc(bp128_M, (uint64_t) sizeof(bp128_t*) * (max_text_len + 2ull + bp128_n_words)))
		bp128_M = t;
	else
		exit(1);

	if (auto t = (bp128_t*)alloc_aligned(bp128_raw_ptr_M, (uint64_t)(this->max_text_len + 2ull + bp128_n_words) * bp128_n_words * sizeof(bp128_t), sizeof(simd128_t)))
		bp128_raw_M = t;
	else
		exit(1);

	for (uint32_t i = 0; i < this->max_text_len + 2; ++i)
		bp128_M[i] = &bp128_raw_M[(uint64_t) i * bp128_n_words];

	for (uint32_t i = 0; i < bp128_n_words; ++i)
	{
		bp128_M[0][i].VP = ~(0ull);
		bp128_M[0][i].VN = 0;

		bp128_M[0][i].HN = 0;
		bp128_M[0][i].HP = 0;
		bp128_M[0][i].D0 = ~(0ull);
	}

	genome_prefetch = (uchar_t*)realloc(genome_prefetch, (uint64_t) sizeof(uchar_t) * (std::max(max_query_len, this->max_text_len) + 8ull)); // fixme: why + 5?
}

// ************************************************************************************
template <instruction_set_t instruction_set>
bool LevMyers128<instruction_set>::dynamicProgramming(
	ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance)
{
	uint32_t j;
	uint32_t seq_len_m1 = seq_len - 1;
	uint32_t red_m = (seq_len_m1) % 128;
	uint64_t red_m_mask = 1ull << (red_m % 64);
	uint64_t red_m_word = red_m / 64;

	uint32_t word_num = seq_len_m1 / 64;
	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;// - bp128_n_words + 1;

	//int middle_point_in_ref = (max_distance_in_ref + seq_len) / 2;
	//int best_dist_from_middle = max_distance_in_ref + seq_len;		// large value

	if (max_distance_in_ref == 0) {
		max_distance_in_ref = max_text_len;
		throw std::runtime_error("max_distance_in_ref == 0");
	}

	if (max_distance_in_ref > max_text_len) {
		throw std::runtime_error("reallocation needed: max_query_len: " + to_string(max_text_len) + " max_dist: " + to_string(max_distance_in_ref));
		reallocBuffers(max_query_len, max_distance_in_ref, 2);
	}

	bp128_t * __restrict curr_bp_M = bp128_M[1];

	// Assumption: seq_len <= 128
	// Genome prefetch
	uchar_t * __restrict ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	//uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;

	*ptr++ = 7;			// No symbol
	PrefetchDecompressed128<instruction_set>(ptr, genome_ptr, max_distance_in_ref + 3, ref_pos & 1);
	ptr = genome_prefetch;

	simd128_t HN = bp128_M[0][0].HN;
	simd128_t HP = bp128_M[0][0].HP;
	simd128_t VN = bp128_M[0][0].VN;
	simd128_t VP = bp128_M[0][0].VP;
	simd128_t D0;

	uint32_t ed_threshold = max_distance_in_ref + 1 - 1 + max_mate_edit_distance; // j = 1

	for (j = 1; j <= max_distance_in_ref + 1; ++j, --ed_threshold)
	{
		uint64_t X1 = bp_PM[*ptr][1];
		++ptr;
		uint64_t X0 = bp_PM[*ptr][0];
		simd128_t X(X0, X1);

		simd128_t permshift_HN = permute2<-1, 0>(HN);
		simd128_t permshift_HP = permute2<-1, 0>(HP);
		permshift_HN >>= 63;
		permshift_HP >>= 63;

		X |= permshift_HN;
		//D0 = (((X & VP) + VP) ^ VP) | X | VN;
		D0 = X & VP;
		D0 += VP;
		D0 ^= VP;
		D0 |= X;
		D0 |= VN;

		//HP = VN | ~(D0 | VP);
		HP = D0 | VP;
		HP = ~HP;
		HP |= VN;
		HN = D0 & VP;

		X = HP << 1;
		X |= permshift_HP;
		VN = D0 & X;

		//VP = (HN << 1) | ~(D0 | X);
		X |= D0;
		X = ~X;
		VP = HN << 1;
		VP |= X;

		VP |= permshift_HN;

		// store values in array
		curr_bp_M->D0 = D0;
		curr_bp_M->HN = HN;
		curr_bp_M->HP = HP;
		curr_bp_M->VN = VN;
		curr_bp_M->VP = VP;

		curr_bp_M++;

//		if (HN.get_bit(red_m))
		if (HN[red_m_word] & red_m_mask)
		{
			if (--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j - word_num;
				//best_dist_from_middle = abs((int) j - middle_point_in_ref);
			}
			else if (curr_ed == min_ed)
			{
/*				if (D0.get_bit(red_m))
				{
					int dist_from_middle = abs((int)j - middle_point_in_ref);
					if (dist_from_middle < best_dist_from_middle)
					{
						best_dist_from_middle = dist_from_middle;
						min_ed_pos = j - word_num;
					}
				}
	*/			
				if (min_ed_pos + 1 == j - word_num) // prefer SNP over INDEL
					min_ed_pos = j - word_num;
			}
		}
		else
		{
//			if (HP.get_bit(red_m))
			if (HP[red_m_word] & red_m_mask)
				curr_ed++;
			else if (curr_ed == min_ed)
			{
				if (min_ed_pos + 1 == j - word_num)
					min_ed_pos = j - word_num;
/*				if (D0.get_bit(red_m))
				{
					int dist_from_middle = abs((int)j - middle_point_in_ref);
					if (dist_from_middle < best_dist_from_middle)
					{
						best_dist_from_middle = dist_from_middle;
						min_ed_pos = j - word_num;
					}
				}*/
			}

			if (curr_ed > ed_threshold)
				break;
		}
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
ref_pos_t LevMyers128<instruction_set>::getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
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
	num_events = 0;
	bool inIndel = false;
	ext_cigar[0] = 0;
	string_reader_t quality_reader(quality, seq_len, dir);

	while (read_pos&& pos) {
		uint32_t word_no = (read_pos - 1) / 64;
		uint64_t word_mask = 1ull << ((read_pos - 1) % 64);

		bp128_t word128 = *bp128_M[pos + word_no];

		d_val = h_val = v_val = c_val;

		if (!(word128.D0.extract(word_no) & word_mask)) { --d_val; }
		if (word128.HP.extract(word_no) & word_mask) { --h_val; }
		if (word128.HN.extract(word_no) & word_mask) { ++h_val; }
		if (word128.VP.extract(word_no) & word_mask) { --v_val; }
		if (word128.VN.extract(word_no) & word_mask) { ++v_val; }

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
	}

	ext_cigar[cigar_pos] = '\0';

	reverse(ext_cigar, ext_cigar + cigar_pos);
	return pos;
}

// ************************************************************************************
#ifdef _DEVELOPMENT_MODE
template <instruction_set_t instruction_set>
void LevMyers128<instruction_set>::save(const std::string& filename)
{
	std::ofstream file(filename);

	for (int i = 0; i < (int)max_text_len; ++i) {
		file << std::dec << i << "," << code2symbol[genome_prefetch[i]];

		file << endl << "D0: " << printVec2uq(bp128_M[i][0].D0) << ","
			<< endl << "HN: " << printVec2uq(bp128_M[i][0].HN) << ","
			<< endl << "HP: " << printVec2uq(bp128_M[i][0].HP) << ","
			<< endl << "VN: " << printVec2uq(bp128_M[i][0].VN) << ","
			<< endl << "VP: " << printVec2uq(bp128_M[i][0].VP) << ","
			<< endl << endl;
	}
}
#endif

// EOF


