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
#include "simd_utils.h"

// ************************************************************************************
// 256-bit general implementation
template <instruction_set_t instruction_set>
size_t CountEOLs256(uchar_t* ptr, size_t size)
{
	if (size < 128)
		return CountEOLs(ptr, size);

	size_t n_eols = 0;

	// Allign the starting pointer
	while (reinterpret_cast<uint64_t>(ptr) % 64)
	{
		n_eols += *ptr++ == 0xA;
		--size;
	}

	// Count in 64-byte packs
	size_t n_packs = size / 64;
	uchar_t* pack_end_ptr = ptr + n_packs * 64;

	Vec32uc v_EOL(0xA);
	Vec32uc v_data;
	Vec32cb v_cmp;

	while (ptr != pack_end_ptr)
	{
		v_data.load_a(ptr);
		v_cmp = v_data == v_EOL;
		n_eols += horizontal_count(v_cmp);
		ptr += 32;

		v_data.load_a(ptr);
		v_cmp = v_data == v_EOL;
		n_eols += horizontal_count(v_cmp);
		ptr += 32;
	}

	// Add remaining 
	uchar_t* end_ptr = ptr + (size - n_packs * 64);

	while (ptr != end_ptr)
		n_eols += *ptr++ == 0xA;

	return n_eols;
}

// ************************************************************************************
// 256-bit general implementation
template <instruction_set_t instruction_set>
size_t CountMismatches256(uchar_t* p, uchar_t* q, size_t size)
{
	Vec32uc v_p, v_q;
	size_t rest = size % 32;
	size_t r = 0;

	r += CountMismatches128<instruction_set>(p, q, rest);
	p + rest;
	q += rest;

	for (size -= rest; size; p += 32, q += 32, size -= 32)
	{
		v_p.load(p);
		v_q.load(q);

		r += horizontal_count(v_p != v_q);
	}

	return r;
}

// ************************************************************************************
// 256-bit general implementation
template <instruction_set_t instruction_set>
void PrefetchDecompressed256(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte)
{
	if (first_single_byte)
	{
		*dest++ = *src++ & 0x0f;
		--size;
	}

	// For simplicity of implementation 1 more byte could be prefetched
	if (size & 1)
		++size;

	size_t rest = size % 32;

	if (rest)
	{
		PrefetchDecompressed64(dest, src, rest, false);
		size -= rest;
		src += rest / 2;
		dest += rest;
	}

	Vec16uc v_in;
	Vec16us v_tmp, v_out;
	Vec16us v_mask(0xf0f);

	for (; size; size -= 32)
	{
		v_in.load(src);
		v_tmp = extend(v_in);
		v_out = (v_tmp << 8) ^ (v_tmp >> 4);
		v_out &= v_mask;

		v_out.store(dest);

		dest += 32;
		src += 16;
	}
}

// EOF

