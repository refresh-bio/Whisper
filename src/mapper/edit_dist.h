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


#ifndef _EDIT_DIST_H
#define _EDIT_DIST_H

#include "../common/defs.h"
#include "../libs/vectorclass.h"

class CEditDist 
{
	uint32_t max_query_len;
	uint32_t max_text_len;
	uint32_t max_ed;

	// Structures for LevDiag algorithm
	uint32_t* prev_dp_ptr;
	uint32_t* curr_dp_ptr;
	uchar_t* genome_prefetch;
	uint32_t cur_offset;

	uchar_t *ref_ptr;
	uint32_t ref_size;

	// Structures for Myers's bit-parallel algorithm (64-bit version)
#ifdef LEV_MYERS_ASSERTIONS
	uint32_t bp_n_words;
	uint64_t **bp_PM;
	uint64_t *raw_bp_PM;

	typedef struct {
		uint64_t D0, VP, VN, HN, HP;} bp_t;
	bp_t *bp_raw_M;
	bp_t **bp_M;

	// Structures for Myers's bit-parallel algorithm (128-bit SSE2 version)
	uint32_t bp128_n_words;
	Vec2uq **bp128_PM;
	Vec2uq *bp128_raw_PM;
	uint32_t bp128_n_chars;
	Vec16c *bp128_chars;

	typedef struct {
		Vec2uq D0, VP, VN, HN, HP;} bp128_t;
	bp128_t *bp128_raw_M;
	bp128_t **bp128_M;

	uint32_t *bp_scores;
#endif

	uint32_t alloc_text_M;

	uint32_t rev_comp_code[8];
	uint32_t raw_code[128];
	uint32_t raw_rev_comp_code[128];
	uint32_t seq_len;

	void allocate();
	void release();

public:
	CEditDist(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed);
	~CEditDist();

	void SetReference(uchar_t *_ref_ptr, uint32_t _ref_size, uint32_t _cur_offset);

	uint32_t LevDiag(uchar_t* pattern_ptr, uint32_t pattern_len, uint32_t position_in_text, 
								uint32_t fixed_left, uint32_t fixed_right, uint32_t max_mismatches, bool dir_gen_flag);

#ifdef LEV_MYERS_ASSERTIONS
	bool LevMyers_PP(uchar_t *seq, uint32_t seq_len, genome_t orientation);
	bool LevMyers_PP_rawseq(uchar_t *seq, uint32_t seq_len, genome_t orientation);
	bool LevMyers(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance);

	bool LevMyers128_PP(uchar_t *seq, uint32_t seq_len, genome_t orientation);
	bool LevMyers128(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance);
#endif
};

#endif

// EOF
