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


#include "simd_utils.h"

// ************************************************************************************
// Basic implementations without SIMD extensions
// ************************************************************************************

// ************************************************************************************
size_t CountEOLs(uchar_t *ptr, size_t size)
{
	size_t n_eols = 0;
	uchar_t *end_ptr = ptr + size;

	switch (size % 8)
	{
	case 8:	n_eols += *ptr++ == 0xA;
	case 7:	n_eols += *ptr++ == 0xA;
	case 6:	n_eols += *ptr++ == 0xA;
	case 5:	n_eols += *ptr++ == 0xA;
	case 4:	n_eols += *ptr++ == 0xA;
	case 3:	n_eols += *ptr++ == 0xA;
	case 2:	n_eols += *ptr++ == 0xA;
	case 1:	n_eols += *ptr++ == 0xA;
	}

	while (ptr != end_ptr)
	{
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
		n_eols += *ptr++ == 0xA;
	}

	return n_eols;
}

// ************************************************************************************
size_t CountMismatches(uchar_t* p, uchar_t* q, size_t size)
{
	size_t r = 0;
	uchar_t *p_end = p + size;

	while (p != p_end)
		r += *p++ != *q++;

	return r;
}

// ************************************************************************************
void PrefetchDecompressed(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte)
{
	if (first_single_byte)
	{
		*dest++ = *src++ &0x0f;
		--size;
	}

	// For simplicity of implementation 1 more byte could be prefetched
	if (size & 1)
		++size;

	for (; size; size -= 2)
	{
		*dest++ = *src >> 4;
		*dest++ = *src++ & 0x0f;
	}
}

// EOF
