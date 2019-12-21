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
size_t CountEOLs64(uchar_t *ptr, size_t size)
{
	return CountEOLs(ptr, size);
}

// ************************************************************************************
size_t CountMismatches64(uchar_t* p, uchar_t* q, size_t size)
{
	return CountMismatches64(p, q, size);
}

// ************************************************************************************
void PrefetchDecompressed64(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte)
{
	PrefetchDecompressed(dest, src, size, first_single_byte);
}

// EOF
