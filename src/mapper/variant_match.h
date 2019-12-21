#pragma once
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

#ifdef ENABLE_VCF_VARIANTS
#include "../common/defs.h"
#include <vector>

using namespace std;

class CVariantMatch 
{
	const uint32_t max_seed_len = 16;
	uchar_t sym2dna[256];
	uchar_t dna2sym[256];

public: 
	CVariantMatch();
	~CVariantMatch();

	bool find_pattern(uchar_t *read, uint32_t read_len, uchar_t *pattern, uint32_t pattern_len, vector<uint32_t> &v_pos);
};

#endif

// EOF