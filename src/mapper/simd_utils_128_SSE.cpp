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


#include "simd_utils_128.hpp"

// ************************************************************************************
// 128-bit SSE specialization
template size_t CountEOLs128<instruction_set_t::sse2>(uchar_t* ptr, size_t size);


// ************************************************************************************
// 128-bit SSE specialization
template size_t CountMismatches128<instruction_set_t::sse2>(uchar_t* p, uchar_t* q, size_t size);

// ************************************************************************************
// 128-bit SSE specialization
template void PrefetchDecompressed128<instruction_set_t::sse2>(uchar_t* dest, uchar_t* src, uint32_t size, bool first_single_byte);

// EOF
