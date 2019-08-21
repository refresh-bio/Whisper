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


#pragma once
#include "LevMyers.h"

// ************************************************************************************
template <instruction_set_t instruction_set>
class LevMyers256 : public LevMyers64 
{
protected:
	// Structures for Myers's bit-parallel algorithm (256-bit AVX2 version)
	typedef Vec4uq simd256_t;

	typedef struct {
		simd256_t D0, VP, VN, HN, HP;
	} bp256_t;

	uint32_t bp256_n_words;
	uchar_t* genome_prefetch;
	void *bp256_raw_ptr_M;
	bp256_t *bp256_raw_M;
	bp256_t **bp256_M;

public:
	LevMyers256(uint32_t _max_text_len, uint32_t _max_ed);
	
	~LevMyers256() {
		free(bp256_M);
		free(bp256_raw_ptr_M);
		free(genome_prefetch);
	}

	virtual bool preprocess(uchar_t *seq, uint32_t seq_len, genome_t orientation) override {
		return LevMyers64::preprocess(seq, seq_len, orientation);
	}

	virtual bool preprocessRawSeq(uchar_t *seq, uint32_t seq_len, genome_t orientation) override {
		return LevMyers64::preprocessRawSeq(seq, seq_len, orientation);
	}

	virtual bool dynamicProgramming(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t max_mate_edit_distance, ref_pos_t &pos, uint32_t &edit_distance) override;

	virtual ref_pos_t getExtCigar(uchar_t *ext_cigar, uchar_t* tmp_ref_sequence, uchar_t* tmp_read_sequence, uchar_t* quality,
		ref_pos_t pos, uint32_t edit_dist, uint32_t seq_len, genome_t dir, const scoring_t &scoring, double& affine_score, uint32_t& num_events) override;

#ifdef _DEVELOPMENT_MODE
	virtual void save(const std::string& filename);
#endif

protected:
	virtual void reallocBuffers(uint32_t _max_query_len, uint32_t _max_text_len, int rounding);
};

// EOF 