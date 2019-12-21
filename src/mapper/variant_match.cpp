// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 2.0
// Date    : 2018-10-22
// License : GNU GPL 3
// *******************************************************************************************

#include "variant_match.h"
#include <algorithm>

using namespace std;

#ifdef ENABLE_VCF_VARIANTS
// *******************************************************************************************
CVariantMatch::CVariantMatch()
{
	fill_n(sym2dna, 256, (uchar_t) 'N');
	sym2dna[sym_code_A] = 'A';
	sym2dna[sym_code_C] = 'C';
	sym2dna[sym_code_G] = 'G';
	sym2dna[sym_code_T] = 'T';

	fill_n(dna2sym, 256, (uchar_t) sym_code_N_ref);
	dna2sym['A'] = sym_code_A;
	dna2sym['C'] = sym_code_C;
	dna2sym['G'] = sym_code_G;
	dna2sym['T'] = sym_code_T;
	dna2sym['N'] = sym_code_N_read;
}

// *******************************************************************************************
CVariantMatch::~CVariantMatch()
{

}

// *******************************************************************************************
bool CVariantMatch::find_pattern(uchar_t *read, uint32_t read_len, uchar_t *pattern, uint32_t pattern_len, vector<uint32_t> &v_pos)
{
	v_pos.clear();

	if (pattern_len > read_len)
		return false;

	auto seed_len = min(max_seed_len, pattern_len);

	uint64_t seed = 0ull;
	uint64_t seed_mask = (~0ull) >> (64 - 4 * seed_len);
	uint64_t read_value = 0ull;

	for (uint32_t i = 0; i < seed_len - 1; ++i)
	{
		seed = (seed << 4) + (uint64_t) pattern[i];
		read_value = (read_value << 4) + (uint64_t) read[i];
	}
	seed = (seed << 4) + pattern[seed_len - 1];

	for (int i = seed_len - 1; i < (int) read_len + (int) seed_len - (int) pattern_len; ++i)
	{
		read_value = ((read_value << 4) + (uint64_t) read[i]) & seed_mask;
	
		if (read_value != seed)
			continue;

		uint32_t matching_len = seed_len;
		for (; matching_len < pattern_len; ++matching_len)
			if (pattern[matching_len] != read[i - (seed_len - 1u) + matching_len])
				break;

		if (matching_len == pattern_len)
			v_pos.push_back((uint32_t) (i - (seed_len - 1)));
	}
	
	return !v_pos.empty();
}

#endif

// EOF
