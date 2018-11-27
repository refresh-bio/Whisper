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


#ifndef _Soft_Clipping_
#define _Soft_Clipping_
#include "../common/defs.h"
#include "../common/utils.h"
#include "../common/types.h"
#include "../common/mmgr.h"
#include <vector>

using namespace std;

// ************************************************************************************
class HashTable 
{
public:
	static const uint32_t EMPTY_CELL;
	static const uint32_t hash_len;

	vector<uint32_t> data;
	vector<uint32_t> data_log;

	uint32_t ht_size;
	uint32_t ht_mask;
	uint64_t ht_mult = 0x5555555555555555ull;

	const size_t size() const { return ht_size;  }

	void resize(uint32_t size) {
		data.resize(size, EMPTY_CELL);
		ht_size = size;
		ht_mask = size - 1;
	}

	void reset() {
		for (auto x : data_log)
			data[x] = EMPTY_CELL;
		data_log.clear();
	}

	uint32_t get_pos(uint64_t key)
	{
		return (key * ht_mult) & ht_mask;
	}

	void insert(uint64_t key, uint32_t val) {
		uint32_t pos = get_pos(key);
		while (data[pos] != EMPTY_CELL)
			pos = (pos + 1) & ht_mask;

		data[pos] = val;
		data_log.emplace_back(pos);
	}

	uint32_t atPosition(uint32_t pos) const {
		return data[pos];
	}

	static uint64_t increment(uint64_t key, uchar_t symbol) {
		return ((key << 2) + (symbol)) & ((1ull << 2 * hash_len) - 1);
	}
};

// ************************************************************************************
class CSoftClipping 
{
protected:
	uint32_t max_query_len;
	uint32_t max_text_len;

	uint32_t cur_offset;
	uchar_t *ref_ptr;
	uint32_t ref_size;

	uint32_t rev_comp_code[8];
	uint32_t raw_code[128];
	uint32_t raw_rev_comp_code[128];
	uint32_t alloc_text_M;
	uint32_t seq_len;
	char code2symbol[8];

	vector<uchar_t> query;

	uchar_t* genome_prefetch;
	vector<uint32_t> checked_pos;
	vector<uint32_t> checked_pos_log;
	uint32_t query_len;

	const uint32_t empty_pos_in_query = 1u << 30;
	HashTable ht;

	void clear_ckecked_pos()
	{
		for (auto x : checked_pos_log)
			checked_pos[x] = empty_pos_in_query;
		checked_pos_log.clear();
	}

public:
	CSoftClipping(uint32_t _max_query_len, uint32_t _max_text_len);
	~CSoftClipping();

	void setReference(uchar_t *_ref_ptr, uint32_t _ref_size, uint32_t _cur_offset);
	bool preprocess(uchar_t *seq, uint32_t seq_len, genome_t orientation);
	bool preprocessRawSeq(uchar_t *seq, uint32_t seq_len, genome_t orientation);
	bool match(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t min_exact_len, ref_pos_t &pos, uint32_t &left_clipping, uint32_t &right_clipping,
		uint32_t &left_match, uint32_t &right_match, int32_t &del_size);

	/*
	ExtCigar format description:
	-match: .
	-mismatch: <reference_symbol>
	-insertion: #<read_symbol>
	-deletion: ^<reference symbol>
	-clipping: $
	*/
	void getExtCigar(uchar_t *ext_cigar, const uchar_t* tmp_ref_sequence, const uchar_t* tmp_read_sequence, ref_pos_t pos, uint32_t seq_len, 
		uint32_t left_clipping, uint32_t right_clipping, uint32_t left_match, uint32_t right_match, int32_t del_size,
		const scoring_t &scoring, double& affine_score);

	void clipIllegalPositions(mapping_desc_t& mapping, const scoring_t& scoring, const std::vector<seq_desc_t>& seq_desc, CMemoryPool<uchar_t>& mp_ext_cigar);
	void clipLowQualityPositions(mapping_desc_t& mapping, const uchar_t* quality);
	void clipBoundaryPositions(mapping_desc_t& mapping, const scoring_t& scoring);
	void clipOverlappingPairs(mapping_pair_t& p, const scoring_t& scoring, CMemoryPool<uchar_t>& mp_ext_cigar);

protected:
	void clipLeft(mapping_desc_t& mapping, int readSymbolsCount, const scoring_t& scoring);
	void clipRight(mapping_desc_t& mapping, int readSymbolsCount, const scoring_t& scoring);
};

#endif

// EOF
