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

// EOF
