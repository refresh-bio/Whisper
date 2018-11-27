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

	// Structures for Myers's bit-parallel algorithm (128-bit SSE2 version)
	uint32_t bp128_n_words;
	uchar_t* genome_prefetch;
	void *bp128_raw_ptr_M;
	bp128_t *bp128_raw_M;
	bp128_t **bp128_M;

public:
	LevMyers128(uint32_t _max_query_len, uint32_t _max_text_len, uint32_t _max_ed);

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
	virtual void reallocBuffers(uint32_t _max_query_len, uint32_t _max_text_len, int rounding)
	{
		LevMyers64::reallocBuffers(_max_query_len, _max_text_len, rounding);
		
		bp128_n_words = (this->max_query_len + 128) / 128;
		bp128_M = (bp128_t**)realloc(bp128_M, sizeof(bp128_t*) * (max_text_len + 2 + bp128_n_words));
		bp128_raw_M = (bp128_t*)alloc_aligned(bp128_raw_ptr_M, (this->max_text_len + 2 + bp128_n_words) * bp128_n_words * sizeof(bp128_t), sizeof(simd128_t));

		for (uint32_t i = 0; i < this->max_text_len + 2; ++i)
			bp128_M[i] = &bp128_raw_M[i * bp128_n_words];

		for (uint32_t i = 0; i < bp128_n_words; ++i)
		{
			bp128_M[0][i].VP = ~(0ull);
			bp128_M[0][i].VN = 0;

			bp128_M[0][i].HN = 0;
			bp128_M[0][i].HP = 0;
			bp128_M[0][i].D0 = ~(0ull);
		}

		genome_prefetch = (uchar_t*)realloc(genome_prefetch, sizeof(uchar_t) * (std::max(max_query_len, this->max_text_len) + 5)); // fixme: why + 5?
	}
};

// EOF
