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


#include "mapping_core.h"
#include "../libs/asmlib.h"
#include "../common/utils.h"

// ************************************************************************************
// Adjust sizes of groups, according to no. of reads and size of suffix array for the current bin
void CMappingCore::adjust_group_sizes(reads_bin_t bin)
{
	string bin_prefix = prefixes[bin.bin_id];			// bin prefix
	sa_dir_bin_size = (uint32_t)sa_dir->GetSAPartSize(Prefix2Int(bin_prefix), (uint32_t) bin_prefix.size());
	sa_rc_bin_size = (uint32_t)sa_rc->GetSAPartSize(Prefix2Int(bin_prefix), (uint32_t) bin_prefix.size());

	n_reads_in_bin = (uint32_t)bin.count;
	if (!n_reads_in_bin)
		n_reads_in_bin = 1;

	// How many times the number of suffixes is larger than the number of reads in the current bin
	sa_dir_frac = (double) sa_dir_bin_size / n_reads_in_bin;
	sa_rc_frac = (double) sa_rc_bin_size / n_reads_in_bin;
	sa_frac = (sa_dir_frac + sa_rc_frac) / 2;

	sa_dir_skip = (uint32_t)(sa_dir_frac / log(2));
	sa_rc_skip = (uint32_t)(sa_rc_frac / log(2));
	//	sa_dir_skip = sa_dir_frac * 1;
	//	sa_rc_skip  = sa_rc_frac  * 1;

	if (sa_dir_skip < 1)
		sa_dir_skip = 1;
	if (sa_rc_skip < 1)
		sa_rc_skip = 1;

	group_size_substr_length_mapping.clear();
	group_size_substr_length_mapping.push_back(0);					// len 0	
	group_size_substr_length_mapping.push_back(0);					// len 1
	group_size_substr_length_mapping.push_back(0);					// len 2
	group_size_substr_length_mapping.push_back(0);					// len 3	[0, 64)

	group_size_substr_length_mapping.push_back(32);					// len 4	[64, 256)
	group_size_substr_length_mapping.push_back(768);				// len 5	[256, 1024)
	group_size_substr_length_mapping.push_back(4096);				// len 6	[1024, 4096)
	group_size_substr_length_mapping.push_back(12340);				// len 7	[4096, 16384)
																	
	group_size_substr_length_mapping.push_back(~(uint32_t)0);		// len 8 (max_substr_len)

#ifdef USE_128b_SSE_CF_AND_DP
	malicious_group_length = (uint32_t)(no_of_substrings * sa_frac / 2);
#else
	malicious_group_length = no_of_substrings * sa_frac * 1.5;
#endif
}

// ************************************************************************************
// Adjust positions of substrings for fast searching in malicious groups
void CMappingCore::adjust_substr_pos_arr()
{
	substr_pos_all_lengths.clear();

	substr_pos_in_text_all_lengths.clear();
	substr_pos_mapping_all_lengths.clear();

	int32_t mid_pos = (start_pos + end_pos) / 2;		// middle of the segment

	for (int32_t substr_len = 0; substr_len <= (int) max_substr_len; ++substr_len)		// substring length
	{
		substr_pos_all_lengths.push_back(vector<uint32_t>());
		int32_t left = 0;
		int32_t right = min_read_len - substr_len;

		for (uint32_t j = 0; j < no_of_substrings; ++j)
		{
			int32_t lk = left;
			int32_t rk = right;

			if ((int32_t)start_pos - (lk + substr_len) >= rk - (int32_t)end_pos)		// choose longer distance (from left or right side)
			{
				substr_pos_all_lengths[substr_len].push_back(lk);
				left += substr_len;
			}
			else
			{
				substr_pos_all_lengths[substr_len].push_back(rk);
				right -= substr_len;
			}
		}
		reverse(substr_pos_all_lengths[substr_len].begin(), substr_pos_all_lengths[substr_len].end());

		substr_pos_in_text_all_lengths.push_back(vector<int32_t>());
		substr_pos_mapping_all_lengths.push_back(vector<pair<uint32_t, uint32_t>>());
		for (uint32_t j = 0, j_left = 0, j_right = 0; j < no_of_substrings; ++j)
		{
			int32_t range;

			if ((int32_t)substr_pos_all_lengths[substr_len][j] < mid_pos)
			{
				range = min(j_left, stage_major);
				j_left++;
			}
			else
			{
				range = min(j_right, stage_major);
				j_right++;
			}

			for (int32_t k = -range; k <= range; ++k)
			{
				// Check whether the position in text is already present in positions to compute
				int32_t pos_in_text = (int32_t)substr_pos_all_lengths[substr_len][j] + k;
				auto p = find(substr_pos_in_text_all_lengths[substr_len].begin(), substr_pos_in_text_all_lengths[substr_len].end(), pos_in_text);
				if (p == substr_pos_in_text_all_lengths[substr_len].end())
				{
					substr_pos_mapping_all_lengths[substr_len].push_back(make_pair(j, substr_pos_in_text_all_lengths[substr_len].size()));
					substr_pos_in_text_all_lengths[substr_len].push_back(pos_in_text);
				}
				else
				{
					substr_pos_mapping_all_lengths[substr_len].push_back(make_pair(j, p - substr_pos_in_text_all_lengths[substr_len].begin()));
				}
			}
		}
	}
}

// ************************************************************************************
uint32_t CMappingCore::get_substring_value_from_pattern(uint32_t substring_length, uchar_t *pattern, uint32_t substr_pos)
{
	//substr_pos - starting point for evaluatung a substring. Given in a number of characters

	uint32_t substring_value = 0;
	uint32_t i, j;
	uchar_t character;

	//-------------------------------------------------------------------------------------------------
	//			odd position of a character as a starting  point
	//-------------------------------------------------------------------------------------------------
	if (substr_pos & 1)
	{
		for (i = substr_pos / 2, j = 0; j < substring_length / 2; i++, j++)
		{
			character = LO_NIBBLE(pattern[i]);					//the low nibble - the first character

			if (character <= sym_code_T)
				substring_value = (substring_value << 2) + character;
			else
				return (1 << 2 * substring_length);				//return LEN;

			character = HI_NIBBLE(pattern[i + 1]);				//the high nibble - the second character

			if (character <= sym_code_T)
				substring_value = (substring_value << 2) + character;
			else
				return (1 << 2 * substring_length);				//return LEN;
		}

		if (substring_length & 1)
		{
			character = LO_NIBBLE(pattern[i]);					//the low nobble - the last character

			if (character <= sym_code_T)
				substring_value = (substring_value << 2) + character;
			else
				return (1 << 2 * substring_length);				//return LEN;
		}
	}
	//-------------------------------------------------------------------------------------------------
	//			even position of a character as a starting  point
	//-----------------------------------------------------------------------------------------------
	else
	{
		for (i = substr_pos / 2, j = 0; j < substring_length / 2; i++, j++)
		{
			character = HI_NIBBLE(pattern[i]);					//the high nibble - the first character

			if (character <= sym_code_T)
				substring_value = (substring_value << 2) + character;
			else
				return (1 << 2 * substring_length);				//return LEN;

			character = LO_NIBBLE(pattern[i]);					//the low nibble - the second character

			if (character <= sym_code_T)
				substring_value = (substring_value << 2) + character;
			else
				return (1 << 2 * substring_length);				//return LEN;
		}

		if (substring_length & 1)
		{
			character = HI_NIBBLE(pattern[i]);					//the high nibble - the last character

			if (character <= sym_code_T)
				substring_value = (substring_value << 2) + character;
			else
				return (1 << 2 * substring_length);				//return LEN;
		}
	}

	return substring_value;
}

// ************************************************************************************
// Important: substring_length must be <= 8
uint32_t CMappingCore::get_substring_value(uint32_t substring_length, int32_t start_substring_position, uint32_t ref_index_given_in_chars)
{
	//ref_index_given_in_chars = a left-end of a sequence, given in a number of characters
	uint32_t substring_value = 0;
	uint32_t ref_offset = ref_index_given_in_chars + start_substring_position;

	uchar_t *ref_ptr_pos = ref_ptr + ref_offset / 2;

	//-------------------------------------------------------------------------------------------------
	//------------ odd position in the text
	//-------------------------------------------------------------------------------------------------
	if (ref_offset & 1)
	{
		substring_value = dense_symbols[LO_NIBBLE(*ref_ptr_pos++)];
		switch ((substring_length - 1) / 2)
		{
		case 3:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		case 2:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		case 1:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		}

		if (!(substring_length % 2))
			substring_value = (substring_value << 2) + dense_symbols[HI_NIBBLE(*ref_ptr_pos)];

		if (substring_value >> 16)								// if any bit of high part of value is set, some N had to be processed
			return (1 << 2 * substring_length);					//return LEN;
	}
	//-------------------------------------------------------------------------------------------------
	//------------ even position in the text
	//-----------------------------------------------------------------------------------------------
	else
	{
		switch (substring_length / 2)
		{
		case 4:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		case 3:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		case 2:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		case 1:
			substring_value = (substring_value << 4) + dense_symbols[*ref_ptr_pos++];
		}

		if (substring_length % 2)
			substring_value = (substring_value << 2) + dense_symbols[HI_NIBBLE(*ref_ptr_pos)];

		if (substring_value >> 16)								// if any bit of high part of value is set, some N had to be processed
			return (1 << 2 * substring_length);					//return LEN;
	}

	return substring_value;
}

// ************************************************************************************
vector<uint32_t>::iterator CMappingCore::list_merger(vector<pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator>> in_lists, uint32_t min_occ,
	vector<uchar_t> *list_merger_tmp, vector<uint32_t>::iterator out_iter)
{
	uchar_t _min_occ = min_occ;

	sort(in_lists.begin(), in_lists.end(),
		[](const pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator> &x, const pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator> &y)
	{return distance(x.first, x.second) > distance(y.first, y.second); });

	uint64_t sum = 0;
	for (auto p = in_lists.begin(); p != in_lists.end(); ++p)
		sum += p->second - p->first;

	uchar_t *list_merger_data = list_merger_tmp->data();

	for (auto p = in_lists.begin(); p != in_lists.begin() + min_occ - 1; ++p)
		for (auto q = p->first; q != p->second; ++q)
			++((list_merger_data)[*q]);

	for (auto p = in_lists.begin() + min_occ - 1; p != in_lists.end(); ++p)
		for (auto q = p->first; q != p->second; ++q)
			if (++((list_merger_data)[*q]) == _min_occ)
				*out_iter++ = *q;

	// Zeroing
	if (sum * 50 < list_merger_tmp->size())
		for (auto p = in_lists.begin(); p != in_lists.end(); ++p)
			for (auto q = p->first; q != p->second; ++q)
				list_merger_data[*q] = (uchar_t)0;
	else
		A_memset(list_merger_tmp->data(), 0, list_merger_tmp->size());

	return out_iter;
}

// EOF