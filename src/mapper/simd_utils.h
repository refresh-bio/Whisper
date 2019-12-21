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


#ifndef _SIMD_UTILS_H_
#define _SIMD_UTILS_H_

#include "../common/defs.h"
#define MAX_VECTOR_SIZE 256
#include "../common/defs.h"
#pragma warning (disable: 26495 26451 6385)
#include "../libs/vectorclass.h"
#include "../mapper/vector_utils.h"
#pragma warning (default: 26495 26451 6385)

// *** Count EOLS
size_t CountEOLs(uchar_t *ptr, size_t size);
size_t CountEOLs64(uchar_t *ptr, size_t size);

template <instruction_set_t instruction_set> 
	size_t CountEOLs128(uchar_t *ptr, size_t size);
size_t CountEOLs256(uchar_t *ptr, size_t size);

// *** Count Mismatches
size_t CountMismatches(uchar_t* p, uchar_t* q, size_t size);
size_t CountMismatches64(uchar_t* p, uchar_t* q, size_t size);
template <instruction_set_t instruction_set> 
	size_t CountMismatches128(uchar_t* p, uchar_t* q, size_t size);
size_t CountMismatches256(uchar_t* p, uchar_t* q, size_t size);

// *** Prefetch with 2->1 decoding
void PrefetchDecompressed(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte);
void PrefetchDecompressed64(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte);
template <instruction_set_t instruction_set>
	void PrefetchDecompressed128(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte);
void PrefetchDecompressed256(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte);

#endif

// EOF
