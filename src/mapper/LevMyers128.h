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
#include "immintrin.h"

// ************************************************************************************
template <instruction_set_t instruction_set>
class LevMyers128 : public LevMyers64 {
protected:
	typedef Vec2uq simd128_t;

	typedef struct {
		simd128_t D0, VP, VN, HN, HP;
	} bp128_t;

	// Structures for Myers's bit-parallel algorithm (128-bit version)
	uint32_t bp128_n_words;
	uchar_t* genome_prefetch;
	void *bp128_raw_ptr_M;
	bp128_t *bp128_raw_M;
	bp128_t **bp128_M;

public:
	LevMyers128(uint32_t _max_text_len, uint32_t _max_ed);

	~LevMyers128()
	{
		free(bp128_M);
		free(bp128_raw_ptr_M);
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
