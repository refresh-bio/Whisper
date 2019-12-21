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
#pragma warning (disable: 26495 26451 6385)
#include "../libs/vectorclass.h"
#pragma warning (default: 26495 26451 6385)

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
	pair<uint32_t, int32_t> LevDiagVariant(uchar_t* pattern_ptr, uint32_t pattern_len, uint32_t position_in_text,
		uint32_t fixed_left, uint32_t fixed_right, int variant_len, uint32_t max_mismatches, bool dir_gen_flag);
};

#endif

// EOF
