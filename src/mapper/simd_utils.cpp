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

// EOF
