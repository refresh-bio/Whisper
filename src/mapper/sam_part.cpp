// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.1
// Date    : 2018-07-10
// License : GNU GPL 3
// *******************************************************************************************


#include "sam_part.h"
//#include "../libs/libdeflate.h"

// *******************************************************************************************
CSamPart::CSamPart(CParams *params, CObjects *objects)
{
	mapped_part_size = params->sam_part_size;
	gzipped_SAM_level = params->gzipped_SAM_level;

	mapped_part_reserve = 1 << 16;
//	mapped_part_reserve = 2 << 10;

	if (gzipped_SAM_level)		// if gzip will be used, it is necessary to consider the potential data expansion
	{
		// The rule according to compressBound() function from zlib
		mapped_part_reserve += (mapped_part_size >> 12) + (mapped_part_size >> 14) + (mapped_part_size >> 15) + 13 + 20;
	}

	gzip = new CGzipMember(gzipped_SAM_level);

	mp_sam_parts = objects->mp_sam_parts;
}

// *******************************************************************************************
CSamPart::~CSamPart()
{
	delete gzip;
}

// *******************************************************************************************
void CSamPart::Reserve()
{
	mp_sam_parts->Reserve(mapped_part);
	mapped_part_pos = 0;
}

// *******************************************************************************************
bool CSamPart::IsFilled()
{
	return mapped_part_pos + mapped_part_reserve > mapped_part_size;
}

// *******************************************************************************************
sam_block_t CSamPart::GetBlock()
{
	if(gzipped_SAM_level == 0)
		return sam_block_t(mapped_part, mapped_part_pos, sam_results_t::mapped);
	else
	{
#ifndef _DEBUG
		uchar_t *gz_buffer;
		mp_sam_parts->Reserve(gz_buffer);

		size_t comp_size = gzip->Compress(mapped_part, mapped_part_pos, gz_buffer, mapped_part_size);

		auto r = sam_block_t(gz_buffer, comp_size, sam_results_t::mapped);

		mp_sam_parts->Free(mapped_part);

		return r;
#endif
	}
}

// ************************************************************************************
// Append field to line
void CSamPart::AppendPart(uchar_t *s, uint32_t s_len, bool add_tab)
{
	copy_n(s, s_len, mapped_part + mapped_part_pos);
	mapped_part_pos += s_len;

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
// Append field to line
void CSamPart::AppendPart(uchar_t *s, bool add_tab)
{
	while (*s)
		mapped_part[mapped_part_pos++] = *s++;

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
// Append reversed field to line 
void CSamPart::AppendPartRev(uchar_t *s, uint32_t s_len, bool add_tab)
{
	for (int32_t i = s_len - 1; i >= 0; --i)
		mapped_part[mapped_part_pos++] = s[i];

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
// Append reversed and complemented field to line
void CSamPart::AppendPartRevComp(uchar_t *s, uint32_t s_len, bool add_tab)
{
	for (int32_t i = s_len - 1; i >= 0; --i)
	{
		switch (s[i])
		{
		case 'A':
			mapped_part[mapped_part_pos++] = 'T';
			break;
		case 'C':
			mapped_part[mapped_part_pos++] = 'G';
			break;
		case 'G':
			mapped_part[mapped_part_pos++] = 'C';
			break;
		case 'T':
			mapped_part[mapped_part_pos++] = 'A';
			break;
		default:
			mapped_part[mapped_part_pos++] = s[i];
			break;
		}
	}

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
void CSamPart::AppendPart(string s, bool add_tab)
{
	auto s_len = s.length();

	copy_n(s.c_str(), s_len, mapped_part + mapped_part_pos);
	mapped_part_pos += s_len;

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
void CSamPart::AppendPart(int32_t x, bool add_tab)
{
	mapped_part_pos += SInt2PChar(x, mapped_part + mapped_part_pos);

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
void CSamPart::AppendPart(double x, int prec, bool add_tab)
{
	mapped_part_pos += SDouble2PChar(x, prec, mapped_part + mapped_part_pos);

	if (add_tab)
		mapped_part[mapped_part_pos++] = '\t';
}

// ************************************************************************************
// Close the last line
void CSamPart::CloseLine()
{
	mapped_part[mapped_part_pos - 1] = '\n';	// replace last tab by new-line
}


// *******************************************************************************************
//
// *******************************************************************************************
CBamPart::CBamPart(CParams *params, CObjects *objects)
{
	mapped_part_size = params->sam_part_size;
	gzipped_BAM_level = params->gzipped_SAM_level;

	mapped_part_reserve = 2 << 10;

	// The rule according to compressBound() function from zlib
	mapped_part_reserve += (mapped_part_size >> 12) + (mapped_part_size >> 14) + (mapped_part_size >> 15) + 13 + 20;

	bgzf = new CBGZF(gzipped_BAM_level);

	mp_sam_parts = objects->mp_sam_parts;

}

// *******************************************************************************************
CBamPart::~CBamPart()
{
	delete bgzf;
}

// *******************************************************************************************
sam_block_t CBamPart::GetBlock()
{
	uchar_t *bgzf_buffer;
	mp_sam_parts->Reserve(bgzf_buffer);

	size_t comp_size = bgzf->Compress(mapped_part, mapped_part_pos, bgzf_buffer, mapped_part_size);

	auto r = sam_block_t(bgzf_buffer, comp_size, sam_results_t::mapped);

	mp_sam_parts->Free(mapped_part);

	return r;
}

// *******************************************************************************************
void CBamPart::Reserve()
{
	mp_sam_parts->Reserve(mapped_part);
	mapped_part_pos = 0;
}

// *******************************************************************************************
bool CBamPart::IsFilled()
{
	return mapped_part_pos + mapped_part_reserve > mapped_part_size;
}

// *******************************************************************************************
void CBamPart::SetRefId(int32_t _ref_id)
{
	ref_id = _ref_id;
}

// *******************************************************************************************
void CBamPart::SetPos(int32_t _pos)
{
	pos = _pos;
}

// *******************************************************************************************
void CBamPart::SetRefLength(int32_t _ref_length)
{
	ref_length = _ref_length;
}

// *******************************************************************************************
void CBamPart::SetMapq(uint32_t _mapq)
{
	mapq = _mapq;
}

// *******************************************************************************************
void CBamPart::SetFlag(uint32_t _flag)
{
	flag = _flag;
}

// *******************************************************************************************
void CBamPart::SetNextRefId(int32_t _next_ref_id)
{
	next_ref_id = _next_ref_id;
}

// *******************************************************************************************
void CBamPart::SetNextPos(int32_t _next_pos)
{
	next_pos = _next_pos;
}

// *******************************************************************************************
void CBamPart::SetTlen(int32_t _tlen)
{
	tlen = _tlen;
}

// *******************************************************************************************
void CBamPart::SetReadName(uchar_t *_read_name)
{
	read_name = _read_name;

	for (read_name_len = 0; read_name[read_name_len] && read_name[read_name_len] != ' ' && read_name[read_name_len] != '\t'; ++read_name_len)
		;
}

// *******************************************************************************************
void CBamPart::SetCigar(uint32_t *_cigar, uint32_t _cigar_len)
{
	cigar = _cigar;
	cigar_len = _cigar_len;
}

// *******************************************************************************************
void CBamPart::SetSeq(uchar_t *_seq, uint32_t _seq_len)
{
	seq = _seq;
	seq_len = _seq_len;
}

// *******************************************************************************************
void CBamPart::SetQual(uchar_t *_qual)
{
	qual = _qual;
}

// *******************************************************************************************
void CBamPart::SetDirection(genome_t _direction)
{
	direction = _direction;
}

// *******************************************************************************************
void CBamPart::AddAuxInt(array<uchar_t, 2> _tag, int32_t _value)
{
	v_aux_int.push_back(make_pair(_tag, _value));
}

// *******************************************************************************************
void CBamPart::AddAuxString(array<uchar_t, 2> _tag, uchar_t *_value)
{
	v_aux_string.push_back(make_pair(_tag, _value));
}

// *******************************************************************************************
void CBamPart::CloseLine()
{
	uint32_t record_start_pos = mapped_part_pos;

	mapped_part_pos += 4;		// skip space for total length of alignment
	
	StoreInt32LSB(mapped_part + mapped_part_pos, ref_id, 4);	
	mapped_part_pos += 4;

	StoreInt32LSB(mapped_part + mapped_part_pos, pos - 1, 4);
	mapped_part_pos += 4;

	StoreUIntLSB(mapped_part + mapped_part_pos, read_name_len + 1, 1);
	mapped_part_pos++;

	StoreUIntLSB(mapped_part + mapped_part_pos, mapq, 1);
	mapped_part_pos++;

	StoreUIntLSB(mapped_part + mapped_part_pos, reg2bin(pos-1, pos-1+ref_length), 2);
	mapped_part_pos += 2;

	StoreUIntLSB(mapped_part + mapped_part_pos, cigar_len, 2);
	mapped_part_pos += 2;

	StoreUIntLSB(mapped_part + mapped_part_pos, flag, 2);
	mapped_part_pos += 2;

	StoreUIntLSB(mapped_part + mapped_part_pos, seq_len, 4);
	mapped_part_pos += 4;

	StoreInt32LSB(mapped_part + mapped_part_pos, next_ref_id, 4);
	mapped_part_pos += 4;

	StoreInt32LSB(mapped_part + mapped_part_pos, next_pos - 1, 4);
	mapped_part_pos += 4;

	StoreInt32LSB(mapped_part + mapped_part_pos, tlen, 4);
	mapped_part_pos += 4;

	for (int i = 0; i < read_name_len; ++i)
		mapped_part[mapped_part_pos++] = read_name[i];
	mapped_part[mapped_part_pos++] = 0;

	for (uint32_t i = 0; i < cigar_len; ++i)
	{
		StoreUIntLSB(mapped_part + mapped_part_pos, cigar[i], 4);
		mapped_part_pos += 4;
	}

	// Store SEQ
	if (direction == genome_t::direct)
	{
		for (int i = 0; i < seq_len / 2; ++i)
			mapped_part[mapped_part_pos++] = (a_seq_code_dir[seq[2 * i]] << 4) + a_seq_code_dir[seq[2 * i + 1]];
		if (seq_len & 1)
			mapped_part[mapped_part_pos++] = a_seq_code_dir[seq[seq_len - 1]] << 4;
	}
	else
	{
		for (int i = 0; i < seq_len / 2; ++i)
			mapped_part[mapped_part_pos++] = (a_seq_code_rc[seq[seq_len - 2 * i - 1]] << 4) + a_seq_code_rc[seq[seq_len - (2 * i + 1) - 1]];
		if (seq_len & 1)
			mapped_part[mapped_part_pos++] = a_seq_code_rc[seq[0]] << 4;
	}

	if(direction == genome_t::direct)
		for (int i = 0; i < seq_len; ++i)
			mapped_part[mapped_part_pos++] = qual[i] - 33;
	else
		for (int i = seq_len-1; i >= 0; --i)
			mapped_part[mapped_part_pos++] = qual[i] - 33;

	for (auto x : v_aux_int)
	{
		mapped_part[mapped_part_pos++] = x.first[0];
		mapped_part[mapped_part_pos++] = x.first[1];
		mapped_part[mapped_part_pos++] = 'i';
		StoreInt32LSB(mapped_part + mapped_part_pos, x.second, 4);
		mapped_part_pos += 4;
	}

	for (auto &x : v_aux_string)
	{
		mapped_part[mapped_part_pos++] = x.first[0];
		mapped_part[mapped_part_pos++] = x.first[1];
		mapped_part[mapped_part_pos++] = 'Z';
		for (int i = 0; x.second[i]; ++i)
			mapped_part[mapped_part_pos++] = x.second[i];
		mapped_part[mapped_part_pos++] = 0;
	}

	StoreUIntLSB(mapped_part + record_start_pos, mapped_part_pos - record_start_pos - 4, 4);

	v_aux_int.clear();
	v_aux_string.clear();
}

// *******************************************************************************************
// Function copied from SAM documentation (22 May 2018 version)
int CBamPart::reg2bin(int beg, int end)
{
	--end;

	if (beg >> 14 == end >> 14) 
		return ((1 << 15) - 1) / 7 + (beg >> 14);
	if (beg >> 17 == end >> 17) 
		return ((1 << 12) - 1) / 7 + (beg >> 17);
	if (beg >> 20 == end >> 20) 
		return ((1 << 9) - 1) / 7 + (beg >> 20);
	if (beg >> 23 == end >> 23) 
		return ((1 << 6) - 1) / 7 + (beg >> 23);
	if (beg >> 26 == end >> 26) 
		return ((1 << 3) - 1) / 7 + (beg >> 26);
	
	return 0;
}

// EOF