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

#ifndef _SAM_PART_H
#define _SAM_PART_H

#include "../common/defs.h"
#include "params.h"

#include <array>

#define STORE_EXTRA_SAM_FIELDS

// *******************************************************************************************
class CSamPart
{
	bool gzipped_SAM_level;

	uchar_t *mapped_part;
	uint64_t mapped_part_size;
	uint64_t mapped_part_pos;
	uint64_t mapped_part_reserve;

	CMemoryPool<uchar_t> *mp_sam_parts;

	CGzipMember *gzip;

public:
	CSamPart(CParams *params, CObjects *objects);
	~CSamPart();

	sam_block_t GetBlock();
	void Reserve();
	bool IsFilled();

	void AppendPart(uchar_t *s, uint32_t s_len, bool add_tab = true);
	void AppendPartRev(uchar_t *s, uint32_t s_len, bool add_tab = true);
	void AppendPartRevComp(uchar_t *s, uint32_t s_len, bool add_tab = true);
	void AppendPart(string s, bool add_tab = true);
	void AppendPart(uchar_t* s, bool add_tab = true);
	void AppendPart(int32_t x, bool add_tab = true);
	void AppendPart(double x, int prec, bool add_tab = true);
	void CloseLine();
};

// *******************************************************************************************
class CBamPart
{
	const array<uint8_t, 128> a_seq_code_dir{
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15,  1, 15,  2, 15, 15, 15,  4, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15,  8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
	const array<uint8_t, 128> a_seq_code_rc{ 
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15,  8, 15,  4, 15, 15, 15,  2, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15,  1, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
	
	uint32_t gzipped_BAM_level;

	uchar_t *mapped_part;
	uint64_t mapped_part_size;
	uint64_t mapped_part_pos;
	uint64_t mapped_part_reserve;

	CMemoryPool<uchar_t> *mp_sam_parts;

	CBGZF *bgzf;

	int32_t ref_id;
	int32_t pos;
	int32_t ref_length;
	uint32_t mapq;
	uint32_t flag;
	int32_t next_ref_id;
	int32_t next_pos;
	int32_t tlen;
	uchar_t *read_name;
	uint32_t read_name_len;
	uint32_t *cigar;
	uint32_t cigar_len;
	uchar_t *seq;
	uint32_t seq_len;
	uchar_t *qual;
	genome_t direction;

	vector<pair<array<uchar_t, 2>, int32_t>> v_aux_int;
	vector<pair<array<uchar_t, 2>, float>> v_aux_float;
	vector<pair<array<uchar_t, 2>, uchar_t*>> v_aux_cstring;
	vector<pair<array<uchar_t, 2>, string>> v_aux_string;

	int reg2bin(int beg, int end);

public:
	CBamPart(CParams *params, CObjects *objects);
	~CBamPart();

	sam_block_t GetBlock();
	void Reserve();
	bool IsFilled();

	void SetRefId(int32_t _ref_id);
	void SetPos(int32_t _pos);
	void SetRefLength(int32_t _ref_length);
	void SetMapq(uint32_t _mapq);
	void SetFlag(uint32_t _flag);
	void SetNextRefId(int32_t _next_ref_id);
	void SetNextPos(int32_t _next_pos);
	void SetTlen(int32_t _tlen);
	void SetReadName(uchar_t *_read_name);
	void SetCigar(uint32_t *_cigar, uint32_t _cigar_len);
	void SetSeq(uchar_t *_seq, uint32_t _seq_len);
	void SetQual(uchar_t *_qual);
	void SetDirection(genome_t direction);

	void AddAuxInt(array<uchar_t, 2> _tag, int32_t _value);
	void AddAuxString(array<uchar_t, 2> _tag, uchar_t *_value);
	void AddAuxString(array<uchar_t, 2> _tag, string _value);
	void AddAuxFloat(array<uchar_t, 2> _tag, float _value);

	void CloseLine();
};

#endif

// EOF