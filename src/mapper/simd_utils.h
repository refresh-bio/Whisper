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


#ifndef _SIMD_UTILS_H_
#define _SIMD_UTILS_H_

#include "../common/defs.h"
#define MAX_VECTOR_SIZE 256
#include "../common/defs.h"
#include "../libs/vectorclass.h"
#include "../mapper/vector_utils.h"

size_t CountEOLs(uchar_t *ptr, size_t size);
size_t CountEOLs64(uchar_t *ptr, size_t size);

template <instruction_set_t instruction_set>
size_t CountEOLs128(uchar_t *ptr, size_t size);

size_t CountEOLs256(uchar_t *ptr, size_t size);

#endif

// EOF
