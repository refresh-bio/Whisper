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


#ifndef _LevMyers_
#define _LevMyers_

#include "../common/defs.h"
#include "../common/utils.h"
#include "simd_utils.h"
#include "../mapper/vector_utils.h"
#include <algorithm>

class scoring_t;

// ************************************************************************************
// LevMyers
// ************************************************************************************
class LevMyers {
protected:
	uint32_t max_query_len;
	uint32_t max_text_len;

	uint32_t cur_offset;
	uchar_t* ref_ptr;
	uint32_t ref_size;

	uint32_t rev_comp_code[16];
	uint32_t raw_code[128];
	uint32_t raw_rev_comp_code[128];
	uint32_t alloc_text_M;
	uint32_t seq_len;
	uchar_t last_symbol;
	char code2symbol[16];
	   
public:
	LevMyers(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed);
	virtual ~LevMyers();

	virtual void setReference(uchar_t* _ref_ptr, uint32_t _ref_size, uint32_t _cur_offset);
	virtual bool preprocess(uchar_t* seq, uint32_t seq_len, genome_t orientation) = 0;
	virtual bool preprocessRawSeq(uchar_t* seq, uint32_t seq_len, genome_t orientation) = 0;
	virtual bool dynamicProgramming(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t& pos, uint32_t& edit_distance) = 0;

	enum class matrix_dir_t { H, V, D };

	/*
	ExtCigar format description:
	-match: .
	-mismatch: <reference_symbol>
	-insertion: #<read_symbol>
	-deletion: ^<reference symbol>
	*/
	virtual ref_pos_t getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality, 
		ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) = 0;

#ifdef _DEVELOPMENT_MODE
	virtual void save(const std::string& filename) = 0;
#endif
};

// ************************************************************************************
// LevMyers64
// ************************************************************************************
class LevMyers64 : public LevMyers {
protected:
	typedef struct {
		uint64_t D0, VP, VN, HN, HP;
	} bp_t;
	
	// Structures for Myers's bit-parallel algorithm (64-bit version)
	uint32_t bp_n_words;
	uint64_t **bp_PM;
	uint64_t *raw_bp_PM;

	void *bp_raw_ptr_M;
	bp_t *bp_raw_M;
	bp_t **bp_M;

public:
	LevMyers64(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed);
	LevMyers64(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed, int rounding);
	virtual ~LevMyers64();
	
	virtual bool preprocess(uchar_t *seq, uint32_t seq_len, genome_t orientation) override;
	virtual bool preprocessRawSeq(uchar_t *seq, uint32_t seq_len, genome_t orientation) override;
	virtual bool dynamicProgramming(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance) override;
	virtual ref_pos_t getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
		ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) override;

#ifdef _DEVELOPMENT_MODE
	virtual void save(const std::string& filename);
#endif

protected:
	virtual void reallocBuffers(uint32_t _max_query_len, uint32_t _max_text_len, int rounding);
};

#endif

// EOF
