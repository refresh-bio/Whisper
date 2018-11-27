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
#include "LevMyers128Native.h"
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
template<>
LevMyers128<instruction_set_t::avx>::LevMyers128(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed)
	: 
	LevMyers64(_max_query_len, _max_text_len, _max_ed, 2),
	genome_prefetch(nullptr),
	bp128_raw_ptr_M(nullptr),
	bp128_M(nullptr)
{
	reallocBuffers(_max_query_len, _max_text_len, 2);
}

// ************************************************************************************
template<>
bool LevMyers128<instruction_set_t::avx>::dynamicProgramming(
	ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance) 
{
	uint32_t j;
	uint32_t red_m = (seq_len - 1) % 128;
	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;// - bp128_n_words + 1;

	if (max_distance_in_ref == 0) {
		max_distance_in_ref = max_text_len;
		throw std::runtime_error("max_distance_in_ref == 0");
	}

	if (max_distance_in_ref > max_text_len) {
		reallocBuffers(max_query_len, max_distance_in_ref, 2);
		throw std::runtime_error("reallocation needed");
	}

	bp128_t * __restrict curr_bp_M = bp128_M[1];

	// Assumption: seq_len <= 128
	// Genome prefetch
	uchar_t * __restrict ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;

	*ptr++ = 7;			// No symbol
	if (ref_pos & 1)
		*ptr++ = *genome_ptr++ & 0x0f;
	for (uint32_t i = 0; i < text_len_div2; ++i)
	{
		*ptr++ = *genome_ptr >> 4;
		*ptr++ = *genome_ptr++ & 0x0f;
	}

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
		simd128_t X (X0, X1);
		
		simd128_t permshift_HN = permute2uq<-1, 0>(HN);
		simd128_t permshift_HP = permute2uq<-1, 0>(HP);
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

		if (HN.get_bit(red_m))
		{
			if (--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j - (seq_len / 64); 
			}
			else if (curr_ed == min_ed)
				if (min_ed_pos + 1 == j - (seq_len / 64)) // prefer SNP over INDEL
					min_ed_pos = j - (seq_len / 64); 
		}
		else
		{
			if (HP.get_bit(red_m))
				curr_ed++;
			else if(curr_ed == min_ed)
				if (min_ed_pos + 1 == j - (seq_len / 64))
					min_ed_pos = j - (seq_len / 64); 

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
template<>
ref_pos_t LevMyers128<instruction_set_t::avx>::getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
	ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) 
{
	uint32_t read_pos = seq_len;
	uint32_t cigar_pos = 0;

	const char *decode = "ACGTNNXX";

	uint32_t h_val, v_val, d_val, c_val;

	c_val = edit_dist;
	affine_score = 0;
	num_events = 0;
	bool inIndel = false;
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

			const auto& word128 = *bp128_M[pos + word_no];

			d_val = h_val = v_val = c_val;

			if (!(word128.D0.extract(word_no) & word_mask)) { --d_val; }
			if (word128.HP.extract(word_no) & word_mask) { --h_val; }
			if (word128.HN.extract(word_no) & word_mask) { ++h_val; }
			if (word128.VP.extract(word_no) & word_mask) { --v_val; }
			if (word128.VN.extract(word_no) & word_mask) { ++v_val; }

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
		if (inIndel) {
			affine_score += scoring.gap_extend;
		}
		else {
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
template<>
void LevMyers128<instruction_set_t::avx>::save(const std::string& filename)
{
	std::ofstream file(filename);

	for (int i = 0; i < (int) max_text_len; ++i) {
		file << std::dec << i << "," << code2symbol[genome_prefetch[i]];


		file << endl << "D0: " << printVec2uq(bp128_M[i][0].D0) << ","
			<< endl << "HN: " << printVec2uq(bp128_M[i][0].HN) << ","
			<< endl << "HP: " << printVec2uq(bp128_M[i][0].HP) << ","
			<< endl << "VN: " << printVec2uq(bp128_M[i][0].VN) << ","
			<< endl << "VP: " << printVec2uq(bp128_M[i][0].VP) << ","
		/*	<< endl 
			<< endl << "X0: " << printVec2uq(bp128_M[i][0].X0)
			<< endl << "X1: " << printVec2uq(bp128_M[i][0].X1)
			<< endl << "X2: " << printVec2uq(bp128_M[i][0].X2)
			<< endl << "X3: " << printVec2uq(bp128_M[i][0].X3)
*/
			<< endl << endl;
	}
}
#endif

// ************************************************************************************
template<>
LevMyers128Native<instruction_set_t::avx>::LevMyers128Native(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed)
	:
	LevMyers64(_max_query_len, _max_text_len, _max_ed, 2),
	genome_prefetch(nullptr),
	bp128_raw_ptr_M(nullptr),
	bp128_M(nullptr)
{
	reallocBuffers(_max_query_len, _max_text_len, 2);
}

// ************************************************************************************
template<>
bool LevMyers128Native<instruction_set_t::avx>::dynamicProgramming(
	ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance) 
{
	uint32_t j;
//	uint32_t red_m = (seq_len - 1) % 128;
	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;// - bp128_n_words + 1;

	uint32_t red_m1 = (seq_len - 1) / 8;
	uint32_t red_m2 = 1 << ((seq_len - 1) % 8);

	if (max_distance_in_ref == 0) {
		max_distance_in_ref = max_text_len;
		throw std::runtime_error("max_distance_in_ref == 0");
	}

	if (max_distance_in_ref > max_text_len) {
		reallocBuffers(max_query_len, max_distance_in_ref, 2);
		throw std::runtime_error("reallocation needed");
	}

	bp128_t * __restrict curr_bp_M = bp128_M[1];

	// Assumption: seq_len <= 128
	// Genome prefetch
	uchar_t * __restrict ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;

	*ptr++ = 7;			// No symbol
	if (ref_pos & 1)
		*ptr++ = *genome_ptr++ & 0x0f;
	for (uint32_t i = 0; i < text_len_div2; ++i)
	{
		*ptr++ = *genome_ptr >> 4;
		*ptr++ = *genome_ptr++ & 0x0f;
	}

	ptr = genome_prefetch;

	simd128_t HN = bp128_M[0][0].HN;
	simd128_t HP = bp128_M[0][0].HP;
	simd128_t VN = bp128_M[0][0].VN;
	simd128_t VP = bp128_M[0][0].VP;
	simd128_t D0;

	simd128_t val_63 = _mm_cvtsi64_si128(63);
	simd128_t val_1 = _mm_cvtsi64_si128(1);
	simd128_t val_ones = _mm_set1_epi32(-1);

	simd128_u_t u;

//	uint8_t &mask_u = u.byte[red_m1];

	uint32_t ed_threshold = max_distance_in_ref + 1 - 1 + max_mate_edit_distance; // j = 1

	for (j = 1; j <= max_distance_in_ref + 1; ++j, --ed_threshold)
	{
		uint64_t X1 = bp_PM[*ptr][1];
		++ptr;
		uint64_t X0 = bp_PM[*ptr][0];
//		simd128_t X(X0, X1);
		simd128_t X = _mm_set_epi64x(X1, X0);

//		simd128_t permshift_HN = permute2uq<-1, 0>(HN);
		simd128_t permshift_HN = _mm_slli_si128(HN, 8);
//		simd128_t permshift_HP = permute2uq<-1, 0>(HP);
		simd128_t permshift_HP = _mm_slli_si128(HP, 8);

//		permshift_HN >>= 63;
		permshift_HN = _mm_srl_epi64(permshift_HN, val_63);
//		permshift_HP >>= 63;
		permshift_HP = _mm_srl_epi64(permshift_HP, val_63);

//		X |= permshift_HN;
		X = _mm_or_si128(X, permshift_HN);
		//D0 = (((X & VP) + VP) ^ VP) | X | VN;
//		D0 = X & VP;
		D0 = _mm_and_si128(X, VP);
//		D0 += VP;
		D0 = _mm_add_epi64(D0, VP);
//		D0 ^= VP;
		D0 = _mm_xor_si128(D0, VP);
//		D0 |= X;
		D0 = _mm_or_si128(D0, X);
//		D0 |= VN;
		D0 = _mm_or_si128(D0, VN);

		//HP = VN | ~(D0 | VP);
//		HP = D0 | VP;
		HP = _mm_or_si128(D0, VP);
//		HP = ~HP;
		HP = _mm_xor_si128(HP, val_ones);
//		HP |= VN;
		HP = _mm_or_si128(HP, VN);
//		HN = D0 & VP;
		HN = _mm_and_si128(D0, VP);

//		X = HP << 1;
		X = _mm_sll_epi64(HP, val_1);
//		X |= permshift_HP;
		X = _mm_or_si128(X, permshift_HP);

//		VN = D0 & X;
		VN = _mm_and_si128(D0, X);

		//VP = (HN << 1) | ~(D0 | X);
//		X |= D0;
		X = _mm_or_si128(X, D0);
//		X = ~X;
		X = _mm_xor_si128(X, val_ones);
//		VP = HN << 1;
		VP = _mm_sll_epi64(HN, val_1);
//		VP |= X;
		VP = _mm_or_si128(VP, X);

//		VP |= permshift_HN;
		VP = _mm_or_si128(VP, permshift_HN);

		// store values in array
		curr_bp_M->D0 = D0;
		curr_bp_M->HN = HN;
		curr_bp_M->HP = HP;
		curr_bp_M->VN = VN;
		curr_bp_M->VP = VP;

		curr_bp_M++;

		u.x = HN;

		//		if (HN.get_bit(red_m))
//		if (HN[red_m1] & red_m2)
		if(u.byte[red_m1] & red_m2)
//		if (mask_u & red_m2)
		{
			if (--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j - (seq_len / 64); 
			}
		}
		else
		{
			u.x = HP;
//			if (HP[red_m1] & red_m2)
			if (u.byte[red_m1] & red_m2)
//			if (mask_u & red_m2)
				//			if (HP.get_bit(red_m))
				curr_ed++;
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
template<>
ref_pos_t LevMyers128Native<instruction_set_t::avx>::getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
	ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) 
{
	uint32_t read_pos = seq_len;
	uint32_t cigar_pos = 0;

	const char *decode = "ACGTNNXX";

	uint32_t h_val, v_val, d_val, c_val;

	c_val = edit_dist;
	bool inIndel = false;
	affine_score = 0;
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
//			uint64_t word_mask = 1ull << ((read_pos - 1) % 64);

			uint32_t byte_no = (read_pos - 1) / 8;
			uint64_t byte_mask = 1ull << ((read_pos - 1) % 8);

			const auto& word128 = *bp128_M[pos + word_no];

			simd128_u_t u;

			d_val = h_val = v_val = c_val;

// !!! TODO
/*			if (!(word128.D0.extract(word_no) & word_mask)) { --d_val; }
			if (word128.HP.extract(word_no) & word_mask) { --h_val; }
			if (word128.HN.extract(word_no) & word_mask) { ++h_val; }
			if (word128.VP.extract(word_no) & word_mask) { --v_val; }
			if (word128.VN.extract(word_no) & word_mask) { ++v_val; }
			*/
			u.x = word128.D0;
			if (!(u.byte[byte_no] & byte_mask)) { --d_val; }
			u.x = word128.HP;
			if (u.byte[byte_no] & byte_mask) { --h_val; }
			u.x = word128.HN;
			if (u.byte[byte_no] & byte_mask) { ++h_val; }
			u.x = word128.VP;
			if (u.byte[byte_no] & byte_mask) { --v_val; }
			u.x = word128.VN;
			if (u.byte[byte_no] & byte_mask) { ++v_val; }
			

			if (d_val <= h_val && d_val <= v_val) {
				// mismatch
				affine_score += scoring.mismatch; //quality_reader[read_pos];
				inIndel = false;
				++num_events;

				ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos - 1]];
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

				ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos - 1]];
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
		if (inIndel) {
			affine_score += scoring.gap_extend;
		}
		else {
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
template<>
void LevMyers128Native<instruction_set_t::avx>::save(const std::string& filename)
{
	std::ofstream file(filename);

	for (int i = 0; i < (int) max_text_len; ++i) {
		file << std::dec << i << "," << code2symbol[genome_prefetch[i]];


		file << endl << "D0: " << printVec2uq(bp128_M[i][0].D0) << ","
			<< endl << "HN: " << printVec2uq(bp128_M[i][0].HN) << ","
			<< endl << "HP: " << printVec2uq(bp128_M[i][0].HP) << ","
			<< endl << "VN: " << printVec2uq(bp128_M[i][0].VN) << ","
			<< endl << "VP: " << printVec2uq(bp128_M[i][0].VP) << ","
			/*	<< endl
			<< endl << "X0: " << printVec2uq(bp128_M[i][0].X0)
			<< endl << "X1: " << printVec2uq(bp128_M[i][0].X1)
			<< endl << "X2: " << printVec2uq(bp128_M[i][0].X2)
			<< endl << "X3: " << printVec2uq(bp128_M[i][0].X3)
			*/
			<< endl << endl;
	}
}
#endif

// EOF


