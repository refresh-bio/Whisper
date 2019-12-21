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


#include "LevMyers.h"
#include "../common/types.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#define HI_NIBBLE(x)		((x) >> 4)
#define LO_NIBBLE(x)		((x) & 0x0f)

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) ((i) & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

// ************************************************************************************
LevMyers::LevMyers(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed) :
	cur_offset(0),
	ref_ptr(nullptr),
	ref_size(0),
	alloc_text_M(0),
	seq_len(0),
	last_symbol(0)
{
	max_query_len = _max_query_len;
	max_text_len = _max_text_len;

	rev_comp_code[(int32_t) sym_code_A] = sym_code_T;
	rev_comp_code[(int32_t) sym_code_C] = sym_code_G;
	rev_comp_code[(int32_t) sym_code_G] = sym_code_C;
	rev_comp_code[(int32_t) sym_code_T] = sym_code_A;
	rev_comp_code[4] = 4;
	rev_comp_code[5] = 5;
	rev_comp_code[6] = 6;
	rev_comp_code[7] = 7;
	fill_n(rev_comp_code + 8, 8, 7);

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

	code2symbol[(int32_t) sym_code_A] = 'A';
	code2symbol[(int32_t) sym_code_C] = 'C';
	code2symbol[(int32_t) sym_code_G] = 'G';
	code2symbol[(int32_t) sym_code_T] = 'T';
	code2symbol[(int32_t) sym_code_N_read] = 'N';
	code2symbol[(int32_t) sym_code_N_ref] = 'N';
	fill_n(code2symbol + 8, 8, 'N');
}

// ************************************************************************************
LevMyers::~LevMyers() 
{
}

// ************************************************************************************
void LevMyers::setReference(uchar_t *_ref_ptr, uint32_t _ref_size, uint32_t _cur_offset)
{
	ref_ptr = _ref_ptr;
	ref_size = _ref_size;

	cur_offset = _cur_offset;
}

// ************************************************************************************
LevMyers64::LevMyers64(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed, int rounding)
	:
	LevMyers(_max_query_len, _max_text_len, _max_ed),
	bp_PM(nullptr),
	raw_bp_PM(nullptr),
	bp_raw_M(nullptr),
	bp_M(nullptr)
{
	reallocBuffers(_max_query_len, _max_text_len, rounding);
}

// ************************************************************************************
LevMyers64::LevMyers64(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed)
	: LevMyers64(_max_query_len, _max_text_len, _max_ed, 1) 
{
}

// ************************************************************************************
LevMyers64::~LevMyers64() 
{
	free(raw_bp_PM);
	free(bp_PM);
	free(bp_M);
	free(bp_raw_M);
}

// ************************************************************************************
void LevMyers64::reallocBuffers(uint32_t _max_query_len, uint32_t _max_text_len, int rounding)
{
	this->max_query_len = _max_query_len;
	this->max_text_len = _max_text_len;

	if (this->max_text_len < 10)
		exit(1);

	bp_n_words = (this->max_query_len + 64) / 64;
	bp_n_words = (bp_n_words + rounding - 1) / rounding * rounding;
	bp_n_words = max(1u, bp_n_words);

	if (auto t = (uint64_t*)realloc(raw_bp_PM, (uint64_t) sizeof(uint64_t) * (bp_n_words * 16ull)))		// 15 - largest alphabet code
		raw_bp_PM = t;
	else
		exit(1);

	if (auto t = (uint64_t**)realloc(bp_PM, sizeof(uint64_t*) * 16ull))
		bp_PM = t;
	else
		exit(1);

	for (uint32_t i = 0; i < 16; ++i)
		bp_PM[i] = &raw_bp_PM[(uint64_t)i * bp_n_words];

	if (auto t = (bp_t**)realloc(bp_M, sizeof(bp_t*) * (this->max_text_len + 1ull)))
		bp_M = t;
	else
		exit(1);

	if (auto t = (bp_t*)realloc(bp_raw_M, sizeof(bp_t) * (this->max_text_len + 1ull) * bp_n_words))
		bp_raw_M = t;
	else
		exit(1);

	alloc_text_M = (this->max_text_len + 1) * bp_n_words;

	for (uint32_t i = 0; i < this->max_text_len + 1; ++i)
		bp_M[i] = &bp_raw_M[(uint64_t) i * bp_n_words];

	for (uint32_t i = 0; i < bp_n_words; ++i)
	{
		bp_M[0][i].VP = ~(0ull);
		bp_M[0][i].VN = 0;
	}
}

// ************************************************************************************
bool LevMyers64::preprocess(uchar_t *seq, uint32_t _seq_len, genome_t orientation)
{
	seq_len = _seq_len;

	// Preprocessing
	fill_n(raw_bp_PM, bp_n_words * 16, 0);

	if (orientation == genome_t::direct)
	{
		for (uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = GET_CHAR_FROM_PATTERN(i, seq);
			bp_PM[c][i / 64] |= 1ull << (i % 64);
		}

		last_symbol = GET_CHAR_FROM_PATTERN((seq_len - 1), seq);
	}
	else
	{
		for (uint32_t i = 0; i < seq_len; ++i)
		{
			uint32_t c = rev_comp_code[GET_CHAR_FROM_PATTERN(seq_len - i - 1, seq)];
			bp_PM[c][i / 64] |= 1ull << (i % 64);
		}

		last_symbol = rev_comp_code[GET_CHAR_FROM_PATTERN(0, seq)];
	}

	return true;
}

// ************************************************************************************
bool LevMyers64::preprocessRawSeq(uchar_t *seq, uint32_t _seq_len, genome_t orientation)
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
bool LevMyers64::dynamicProgramming(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance)
{
	uint32_t j, r;

	uint32_t left_end_index_in_gen = ref_pos;

	//	bp_scores[0] = seq_len;

	// Searching
	uint64_t X;
	uint32_t red_m = (seq_len - 1) % 64;
	uint64_t mask_m = 1ull << red_m;
	auto local_bp_n_words = (seq_len - 1 + 64) / 64;
	
	uint32_t min_ed = seq_len;
	uint32_t min_ed_pos = 0;
	uint32_t curr_ed = seq_len;

	int middle_point_in_ref = (max_distance_in_ref + seq_len) / 2;
	int best_dist_from_middle = max_distance_in_ref + seq_len;		// large value

	if (max_distance_in_ref == 0) {
		max_distance_in_ref = max_text_len;
		throw std::runtime_error("max_distance_in_ref == 0");
	}

	if (max_distance_in_ref > max_text_len) {
		throw std::runtime_error("reallocation needed: max_query_len: " + to_string(max_query_len) + " max_dist: " + to_string(max_distance_in_ref));
		reallocBuffers(max_query_len, max_distance_in_ref, 2);
	}

	bp_t *curr_bp_M = bp_M[1];
	bp_t *prev_bp_M = bp_M[0];
	
	for (j = 1; j <= max_distance_in_ref; ++j)
	{
		uint32_t c = GET_CHAR_FROM_GENOME(j - 1);

		if (curr_bp_M - bp_raw_M >= alloc_text_M)
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
		for (r = 1; r < local_bp_n_words; ++r)
		{
			curr_bp_M++;
			prev_bp_M++;

			if (curr_bp_M - bp_raw_M >= alloc_text_M)
			{
				cerr << "To far in bp_M!\n";
				cerr << curr_bp_M - bp_raw_M << "   " << alloc_text_M << "\n";
				exit(1);
			}

			X = bp_PM[c][r];
			
			X |= (curr_bp_M - 1)->HN >> 63;
			
			curr_bp_M->D0 = (((X & prev_bp_M->VP) + prev_bp_M->VP) ^ prev_bp_M->VP) | X | prev_bp_M->VN;
			curr_bp_M->HP = prev_bp_M->VN | ~(curr_bp_M->D0 | prev_bp_M->VP);
			curr_bp_M->HN = curr_bp_M->D0 & prev_bp_M->VP;
			X = curr_bp_M->HP << 1;
			
			X |= (curr_bp_M - 1)->HP >> 63;
	
			curr_bp_M->VP = (curr_bp_M->HN << 1) | ~(curr_bp_M->D0 | X);
			curr_bp_M->VP |= (curr_bp_M - 1)->HN >> 63;
			curr_bp_M->VN = curr_bp_M->D0 & X;
		}

		if (curr_bp_M->HP & mask_m)
			curr_ed++;
		else if (curr_bp_M->HN & mask_m)
		{
			if (--curr_ed < min_ed)
			{
				min_ed = curr_ed;
				min_ed_pos = j;
				best_dist_from_middle = abs((int) j - middle_point_in_ref);
			}
			else if (curr_ed == min_ed)
			{
				if (min_ed_pos + 1 == j)
					min_ed_pos = j;
/*				if (curr_bp_M->D0 & mask_m)
				{
					int dist_from_middle = abs((int)j - middle_point_in_ref);
					if (dist_from_middle < best_dist_from_middle)
					{
						best_dist_from_middle = dist_from_middle;
						min_ed_pos = j;
					}
				}*/
			}
		}
		else if (curr_ed == min_ed)
		{
			if (min_ed_pos + 1 == j)
				min_ed_pos = j;
/*			if (curr_bp_M->D0 & mask_m)
			{
				int dist_from_middle = abs((int)j - middle_point_in_ref);
				if (dist_from_middle < best_dist_from_middle)
				{
					best_dist_from_middle = dist_from_middle;
					min_ed_pos = j;
				}
			}*/
		}

		for (r = local_bp_n_words; r < bp_n_words; ++r)
		{
			curr_bp_M++;
			prev_bp_M++;
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
ref_pos_t LevMyers64::getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
	ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events)
{
	uint32_t read_pos = seq_len;
	uint32_t cigar_pos = 0;

	const char *decode = "ACGTNNXX";

	uint32_t h_val, v_val, d_val, c_val;
	matrix_dir_t prev_dir = matrix_dir_t::D;
	matrix_dir_t curr_dir;

	c_val = edit_dist;

	ext_cigar[0] = 0;
	affine_score = 0;
	num_events = 0;
	bool inIndel = false;

	string_reader_t quality_reader(quality, seq_len, dir);

	while (read_pos && pos) {
		uint32_t word_no = (read_pos - 1) / 64;
		uint64_t word_mask = 1ull << ((read_pos - 1) % 64);

		d_val = h_val = v_val = c_val;

		if (!(bp_M[pos][word_no].D0 & word_mask)) { --d_val; }
		if (bp_M[pos][word_no].HP & word_mask) { --h_val; }
		if (bp_M[pos][word_no].HN & word_mask) { ++h_val; }
		if (bp_M[pos][word_no].VP & word_mask) { --v_val; }
		if (bp_M[pos][word_no].VN & word_mask) { ++v_val; }

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

		if(curr_dir == matrix_dir_t::D)
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
		else if(curr_dir == matrix_dir_t::H)
		{
			// deletion
			if (inIndel) {
				affine_score += scoring.gap_del_extend;
			} else {
				affine_score += scoring.gap_del_open;
				++num_events;
				inIndel = true;
			}

			ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos-1]];
			ext_cigar[cigar_pos++] = '^';
			--pos;
			c_val = h_val;
		}
		else 
		{
			// insertion
			if (inIndel) {
				affine_score += scoring.gap_ins_extend;
			} else {
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
void LevMyers64::save(const std::string& filename)
{
	//ofstream file(filename, ios_base::app);
	std::ofstream file(filename);

	for (uint32_t i = 0; i < max_text_len; ++i) {
		file << std::dec << i;

		file << endl << "D0: ";
		for (uint32_t j = 0; j < bp_n_words; ++j) { file << std::hex << std::setw(16) << bp_M[i][j].D0 << ", "; }
		file << endl << "HN: ";
		for (uint32_t j = 0; j < bp_n_words; ++j) { file << std::hex << std::setw(16) << bp_M[i][j].HN << ", "; }
		file << endl << "HP: ";
		for (uint32_t j = 0; j < bp_n_words; ++j) { file << std::hex << std::setw(16) << bp_M[i][j].HP << ", "; }
		file << endl << "VN: ";
		for (uint32_t j = 0; j < bp_n_words; ++j) { file << std::hex << std::setw(16) << bp_M[i][j].VN << ", "; }
		file << endl << "VP: ";
		for (uint32_t j = 0; j < bp_n_words; ++j) { file << std::hex << std::setw(16) << bp_M[i][j].VP << ", "; }
		file << endl;

		file << endl;
	}
}
#endif

// EOF
