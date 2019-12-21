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


#ifndef _INDELS_H
#define _INDELS_H
#include "../common/defs.h"
#include "../common/utils.h"
#include "../common/types.h"
#include "simd_utils.h"
#include "params.h"
#include <vector>
#include <algorithm>

// ************************************************************************************
#if 0
class CBitChars {
	const uint32_t stack_no_words = 64;
	uint64_t stack_data[64];
	uint64_t* highest_word;
	uint64_t* data;
	uint32_t no_words;
	uint64_t hi_word_mask;
	uint32_t alloc_words;
	bool stack_variant;

	void shift_data()
	{
		uint64_t to_move = 0ull;

		if(stack_variant)
		{
			for (uint32_t i = 0; i < no_words; ++i)
			{
				uint64_t top_bits = stack_data[i] >> 60;
				stack_data[i] <<= 4;
				stack_data[i] += to_move;
				to_move = top_bits;
			}

			*highest_word &= hi_word_mask;
		}
		else
		{
			for (uint32_t i = 0; i < no_words; ++i)
			{
				uint64_t top_bits = data[i] >> 60;
				data[i] <<= 4;
				data[i] += to_move;
				to_move = top_bits;
			}

			*highest_word &= hi_word_mask;
		}
	}

public:
	CBitChars(uint32_t size = 1) {
		data = nullptr;
		stack_variant = true;
		no_words = 1;
		alloc_words = 0;

		Reset(size);
	}

	~CBitChars() {
		if(data)
			delete[] data;
	}

	void Clear()
	{
		if(stack_variant)
			for (uint32_t i = 0; i < no_words; ++i)
				stack_data[i] = 0ull;
		else
			for (uint32_t i = 0; i < no_words; ++i)
				data[i] = 0ull;
	}

	void Reset(uint32_t size)
	{
		hi_word_mask = ~0ull >> (4 * (16 - size % 16));
		no_words = (size + 15) / 16;

		if (no_words < 1)
			exit(1);

		if (no_words > alloc_words)
		{
			if (no_words > stack_no_words)
			{
				stack_variant = false;
				if(data)
					delete[] data;
				alloc_words = no_words;
				data = new uint64_t[no_words];
			}
		}

		if (stack_variant)
			highest_word = stack_data + (no_words - 1);
		else
			highest_word = data + (no_words - 1);

		Clear();
	}

	void Insert(uchar_t x)
	{
		shift_data();

		if (stack_variant)
		{
			if (x < 4)
				stack_data[0] += 1ull << x;
		}
		else
		{
			if (x < 4)
				data[0] += 1ull << x;
		}
	}

	void Insert(uchar_t x, uint32_t pos)
	{
		if (x >= 4)
			return;

		if (stack_variant)
			stack_data[pos / 16] += 1ull << (x + 4 * (pos & 0xf));
		else
			data[pos / 16] += 1ull << (x + 4 * (pos & 0xf));
	}

	uint32_t CountMatches(const CBitChars& x)
	{
		uint32_t r = 0;

		if(stack_variant)
			for (uint32_t i = 0; i < no_words; ++i)
				r += pop_count(stack_data[i] & x.stack_data[i]);
		else
			for (uint32_t i = 0; i < no_words; ++i)
				r += pop_count(data[i] & x.data[i]);

		return r;
	}
};
#endif

// ************************************************************************************
class CBitCharsSingleWord {
	uint64_t data;
	uint64_t word_mask;

public:
	CBitCharsSingleWord(uint32_t size = 1) {
		data = 0ull;

		Reset(size);
	}

	~CBitCharsSingleWord() {
	}

	void Clear()
	{
		data = 0ull;
	}

	void Reset(uint32_t size)
	{
		word_mask = ~0ull >> (4 * (16 - size));

		Clear();
	}

	void Insert(uchar_t x)
	{
		if(x < 4)
			data = ((data << 4) + (1ull << x)) & word_mask;
	}

	void Insert(uchar_t x, uint32_t pos)
	{
		if (x >= 4)
			return;

		data += 1ull << (x + 4 * pos);
	}

	uint32_t CountMatches(const CBitCharsSingleWord& x)
	{
		return pop_count(data & x.data);
	}
};

class CIndelMatching
{
	const int min_side_size = 7;
	const int max_no_side_errors = 1;

//	double gap_open;
//	double gap_extend;
	double gap_ins_open;
	double gap_ins_extend;
	double gap_del_open;
	double gap_del_extend;
	double mismatch_score;
	double match_score;
	double clipping_score;
	double min_clipped_factor;

	uchar_t* ref_ptr;
	uint32_t ref_size;
	uint32_t cur_offset;
	uint32_t query_len;

	uint32_t rev_comp_code[16];
	uint32_t raw_code[128];
	uint32_t raw_rev_comp_code[128];
	char code2symbol[16];

	CBitCharsSingleWord bcsw_ref;
	CBitCharsSingleWord bcsw_query;

	vector<uchar_t> query;
	vector<uchar_t> genome;
	uchar_t* genome_prefetch;

	vector<uint32_t> v_left_candidates;
	vector<uint32_t> v_right_candidates;

	vector<uint32_t> v_left_mismatches;
	vector<uint32_t> v_right_mismatches;

	size_t(*ptr_CountMismatches)(uchar_t*, uchar_t*, size_t);
	void(*ptr_PrefetchDecompressed)(uchar_t* , uchar_t*, uint32_t, bool );

	// Small utils
	uint32_t indel_to_mismatch_score(int32_t size);
	uint32_t rescale_max_no_mismatches(uint32_t tested_size, uint32_t max_no_mismatches)
	{
		return tested_size * max_no_mismatches / query_len;
	}

	bool find_candidates(CBitCharsSingleWord &bcsw_query, vector<uint32_t>& v_cand, uint32_t begin, uint32_t end);
	pair<uint32_t, uint32_t> calc_mismatches_with_indel(uint32_t genome_left_pos, uint32_t genome_right_pos, uint32_t query_left_pos, uint32_t query_right_pos, bool left_margin);

	void append_matching_part(uint32_t mp_len, uchar_t* &ec, uint32_t &i_genome, uint32_t &i_query, 
		uint16_t& err_edit_distance, uint16_t& num_events, double& score, CParams* params);
	void append_indel(int indel_len, uchar_t* &ec, uint32_t& i_genome, uint32_t& i_query, 
		uint16_t& err_edit_distance, uint16_t& num_events, double &score, CParams* params);
	void append_clipping(int clipping_len, uchar_t* &ec, uint32_t& i_genome, uint32_t& i_query, 
		uint16_t& err_edit_distance, uint16_t& num_events, double &score, CParams* params);

	bool get_ext_cigar_indel1(mapping_desc_t& mapping_desc, CParams *params);
	bool get_ext_cigar_indel2(mapping_desc_t& mapping_desc, CParams* params);
	bool get_ext_cigar_clipping_mismatches(mapping_desc_t& mapping_desc, CParams* params);
	bool get_ext_cigar_clipping_indel(mapping_desc_t& mapping_desc, CParams* params);
	bool get_ext_cigar_clipping_clipping(mapping_desc_t& mapping_desc, CParams* params);

public:
	CIndelMatching(CParams *params);
	~CIndelMatching();

	void SetReference(uchar_t* _ref_ptr, uint32_t _ref_size, uint32_t _cur_offset);
	
	void prefetch_genome(uint32_t begin, uint32_t size);

	bool Preprocess(uchar_t* seq, uint32_t seq_len, genome_t orientation);
	bool PreprocessRawSeq(uchar_t* seq, uint32_t seq_len, genome_t orientation);
	bool Match(ref_pos_t ref_pos, uint32_t max_indel_size, uint32_t fixed_segment_pos, uint32_t fixed_segment_size, uint32_t max_no_mismatches,
		candidate_mapping_t &candidate_mapping);
	bool GetExtCigar(ref_pos_t ref_pos, int32_t del_size, uchar_t *ext_cigar, uint32_t &no_mismatches);
	bool GetExtCigarNew(mapping_desc_t &mapping_desc, CParams *params);
};

#endif

// EOF
