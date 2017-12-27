// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.0
// Date    : 2017-12-24
// License : GNU GPL 3
// *******************************************************************************************


#include "simd_utils.h"

// ************************************************************************************
// 128-bit AVX2 implementation
template<>
size_t CountEOLs128<instruction_set_t::avx2>(uchar_t *ptr, size_t size)
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
	uchar_t *pack_end_ptr = ptr + n_packs * 64;

	Vec16uc v_EOL(0xA);
	Vec16uc v_data;
	Vec16cb v_cmp;

	while (ptr != pack_end_ptr)
	{
		v_data.load_a(ptr);
		v_cmp = v_data == v_EOL;
		n_eols += horizontal_count(v_cmp);
		ptr += 16;

		v_data.load_a(ptr);
		v_cmp = v_data == v_EOL;
		n_eols += horizontal_count(v_cmp);
		ptr += 16;

		v_data.load_a(ptr);
		v_cmp = v_data == v_EOL;
		n_eols += horizontal_count(v_cmp);
		ptr += 16;

		v_data.load_a(ptr);
		v_cmp = v_data == v_EOL;
		n_eols += horizontal_count(v_cmp);
		ptr += 16;
	}

	// Add remaining 
	uchar_t *end_ptr = ptr + (size - n_packs * 64);

	while (ptr != end_ptr)
		n_eols += *ptr++ == 0xA;

	return n_eols;
}

// EOF
