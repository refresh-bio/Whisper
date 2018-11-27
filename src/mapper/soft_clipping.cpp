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


#include "soft_clipping.h"
#include <algorithm>

#define HI_NIBBLE(x)		((x) >> 4)
#define LO_NIBBLE(x)		((x) & 0x0f)

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) ((i) & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

const uint32_t HashTable::EMPTY_CELL = ~(0u);
const uint32_t HashTable::hash_len = 20;

// ************************************************************************************
CSoftClipping::CSoftClipping(uint32_t _max_query_len, uint32_t _max_text_len)
{
	max_query_len = _max_query_len;
	max_text_len = _max_text_len;

	uint32_t ht_size;
	for (ht_size = 256; ht_size / max_query_len < 50; ht_size *= 2)
		;
	ht.resize(ht_size);

	rev_comp_code[(int32_t)sym_code_A] = sym_code_T;
	rev_comp_code[(int32_t)sym_code_C] = sym_code_G;
	rev_comp_code[(int32_t)sym_code_G] = sym_code_C;
	rev_comp_code[(int32_t)sym_code_T] = sym_code_A;
	rev_comp_code[4] = 4;
	rev_comp_code[5] = 5;
	rev_comp_code[6] = 6;
	rev_comp_code[7] = 7;

	for (int i = 0; i < 128; ++i)
	{
		raw_code[i] = 6;
		raw_rev_comp_code[i] = 6;
	}

	raw_code['A'] = sym_code_A;
	raw_code['C'] = sym_code_C;
	raw_code['G'] = sym_code_G;
	raw_code['T'] = sym_code_T;
	raw_code['N'] = sym_code_N_read;

	raw_rev_comp_code['A'] = sym_code_T;
	raw_rev_comp_code['C'] = sym_code_G;
	raw_rev_comp_code['G'] = sym_code_C;
	raw_rev_comp_code['T'] = sym_code_A;
	raw_rev_comp_code['N'] = sym_code_N_read;

	code2symbol[(int32_t)sym_code_A] = 'A';
	code2symbol[(int32_t)sym_code_C] = 'C';
	code2symbol[(int32_t)sym_code_G] = 'G';
	code2symbol[(int32_t)sym_code_T] = 'T';
	code2symbol[(int32_t)sym_code_N_read] = 'N';
	code2symbol[(int32_t)sym_code_N_ref] = 'N';

	query.resize(max_query_len);

		// FIXME: + 10
	genome_prefetch = new uchar_t[std::max(max_query_len, max_text_len) + 10];

	checked_pos.resize(std::max(max_query_len, max_text_len) + 10);
}

// ************************************************************************************
CSoftClipping::~CSoftClipping()
{
	delete[] genome_prefetch;
}

// ************************************************************************************
void CSoftClipping::setReference(uchar_t *_ref_ptr, uint32_t _ref_size, uint32_t _cur_offset)
{
	ref_ptr = _ref_ptr;
	ref_size = _ref_size;

	cur_offset = _cur_offset;
}

// ************************************************************************************
bool CSoftClipping::preprocess(uchar_t *seq, uint32_t seq_len, genome_t orientation)
{
	ht.reset();
	// Insert query sequence into HT
	uint64_t ht_val = 0;

	query_len = seq_len;

	if (orientation == genome_t::direct)
	{
		uint32_t i;
		for (i = 0; i < ht.hash_len - 1; ++i)
		{
			uint64_t c = GET_CHAR_FROM_PATTERN(i, seq);
			query[i] = (uchar_t) c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
		}
		for (; i < seq_len; ++i)
		{
			uint64_t c = GET_CHAR_FROM_PATTERN(i, seq);
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
			ht.insert(ht_val, i - (ht.hash_len - 1));
		}
	}
	else
	{
		uint32_t i;
		for (i = 0; i < ht.hash_len - 1; ++i)
		{
			uint64_t c = rev_comp_code[GET_CHAR_FROM_PATTERN(seq_len - i - 1, seq)];
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
		}
		for (; i < seq_len; ++i)
		{
			uint64_t c = rev_comp_code[GET_CHAR_FROM_PATTERN(seq_len - i - 1, seq)];
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
			ht.insert(ht_val, i - (ht.hash_len - 1));
		}
	}

	return true;
}

// ************************************************************************************
bool CSoftClipping::preprocessRawSeq(uchar_t *seq, uint32_t seq_len, genome_t orientation)
{
	ht.reset();

	// Insert query sequence into HT
	uint64_t ht_val = 0;

	query_len = seq_len;

	if (orientation == genome_t::direct)
	{
		uint32_t i;
		for (i = 0; i < ht.hash_len - 1; ++i)
		{
			uint64_t c = raw_code[seq[i]];
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
		}
		for (; i < seq_len; ++i)
		{
			uint64_t c = raw_code[seq[i]];
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
			ht.insert(ht_val, i - (ht.hash_len - 1));
		}
	}
	else
	{
		uint32_t i;
		for (i = 0; i < ht.hash_len - 1; ++i)
		{
			uint64_t c = raw_rev_comp_code[seq[seq_len - i - 1]];
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
		}
		for (; i < seq_len; ++i)
		{
			uint64_t c = raw_rev_comp_code[seq[seq_len - i - 1]];
			query[i] = (uchar_t)c;
			ht_val = HashTable::increment(ht_val, (uchar_t) c);
			ht.insert(ht_val, i - (ht.hash_len - 1));
		}
	}

	return true;
}

// ************************************************************************************
bool CSoftClipping::match(ref_pos_t ref_pos, uint32_t max_distance_in_ref, uint32_t min_exact_len, ref_pos_t &pos, uint32_t &left_clipping, uint32_t &right_clipping,
	uint32_t &left_match, uint32_t &right_match, int32_t &del_size)
{
	uint32_t j;
//	uint32_t left_end_index_in_gen = ref_pos;
	
	// Genome prefetch
	uchar_t * __restrict ptr = genome_prefetch;
	uchar_t *genome_ptr = ref_ptr + (ref_pos >> 1);
	uint32_t text_len_div2 = (max_distance_in_ref + 3) / 2;

	*ptr++ = 7;			// No symbol
	if (ref_pos & 1)
		*ptr++ = *genome_ptr++ & 0x0f;
	for (uint32_t i = 0; i < text_len_div2; ++i)
	{
		*ptr++ = *genome_ptr >> 4;
		*ptr++ = *genome_ptr++ & 0x0f;
	}

	clear_ckecked_pos();

	uint64_t key = 0;

	uint32_t best_text_pos = 0; 
	uint32_t best_query_pos = 0;
	uint32_t best_exact_len = 0;

	int32_t best_score = 0;
	int32_t left_exact_len = 0;
	int32_t left_text_pos = 0;
	uint32_t best_left_len = 0;
	uint32_t best_right_len = 0;

	for (j = 1; j < ht.hash_len; ++j)
		key = HashTable::increment(key, genome_prefetch[j]);

	int ht_hits = 0;
	int ht_first_hits = 0;

	for (; j <= max_distance_in_ref; ++j)
	{
		key = HashTable::increment(key, genome_prefetch[j]);
		uint32_t pos_in_query; // 0-based position in query
		uint32_t ht_pos = ht.get_pos(key);

		bool first = true;

		while ((pos_in_query = ht.atPosition(ht_pos)) != ht.EMPTY_CELL)
		{
			++ht_hits;
			if (first) { ++ht_first_hits;  }
			
			uint32_t pos_in_text = j - (ht.hash_len - 1); // j indicates end of the analyzed fragment
			uint32_t pos_in_text_of_query_front = pos_in_text - pos_in_query; 

			if (pos_in_text_of_query_front <= max_distance_in_ref)
			{
				int32_t len = 0;
				if (checked_pos[pos_in_text_of_query_front] + 1 != pos_in_query)
				{
					uchar_t* qq = query.data() + pos_in_query;
					uchar_t* gg = genome_prefetch + pos_in_text;
					int len_range = query_len - pos_in_query;

					for (; len < len_range && *qq == *gg; ++len, ++qq, ++gg)
						;

					if (len >= (int) min_exact_len)
					{
						if (len > best_score)		// check if the match is the longest so far
						{
							best_score = len;
							best_query_pos = pos_in_query;
							best_text_pos = pos_in_text_of_query_front;
							best_left_len = 0;
							best_right_len = 0;
							best_exact_len = len;
						}

						if (pos_in_query == 0)		// if the match is at the left boundary, store this information
						{
							left_exact_len = len;
							left_text_pos = pos_in_text_of_query_front;
						}
						else if (pos_in_query + len == query_len && left_exact_len > 0)		// the match is at the right boundary and there is a left-boundary match - check long indel
						{
							int32_t text_dist = (int32_t)(pos_in_text + len) - (int32_t)left_text_pos;

							if (left_exact_len + len >= (int32_t) query_len)		// deletion
							{
								left_exact_len -= left_exact_len + len - query_len;

								if (text_dist > (int32_t) query_len)
								{
									int32_t score = (int32_t)query_len - (text_dist - (int32_t)query_len);

									if (score > best_score)
									{
										best_score = score;
										best_left_len = left_exact_len;
										best_right_len = len;
										best_query_pos = 0;
										best_text_pos = left_text_pos;
										left_exact_len = 0;
										best_exact_len = 0;
									}
								}
							}
							else if (left_exact_len + len == text_dist)	// insertion
							{
								int32_t score = text_dist;
								if (score > best_score)
								{
									best_score = score;
									best_left_len = left_exact_len;
									best_right_len = len;
									best_query_pos = 0;
									best_text_pos = left_text_pos;
									left_exact_len = 0;
									best_exact_len = 0;
								}
							}
						}
					}
				}
				
				if (len >= (int32_t) min_exact_len)											// Do not mark for fake matches 
				{
					checked_pos[pos_in_text_of_query_front] = pos_in_query;			// Mark that this position in text was already verified
					checked_pos_log.emplace_back(pos_in_text_of_query_front);
				}
				first = false;
			}

			if (++ht_pos >= ht.size()) { ht_pos = 0; }
		}
	}

	if (best_score < (int32_t) min_exact_len)
		return false;

	if (best_exact_len)		// left - right clipping
	{
		pos = best_text_pos - 1; // fixme: to omit dummy symbol from the beginning?
		left_clipping = best_query_pos;
		right_clipping = query_len - best_query_pos - best_exact_len;
		del_size = 0;
		left_match = 0;
		right_match = 0;
	}
	else					// indel
	{
		pos = best_text_pos - 1;
		left_clipping = 0;
		right_clipping = 0;
		left_match = best_left_len;
		right_match = best_right_len;
		if (best_left_len + best_right_len == query_len)
			del_size = query_len - best_score;
		else
			del_size = -(int32_t)(query_len - left_match - right_match);
	}

	return true;
}

// ************************************************************************************
void CSoftClipping::getExtCigar(
	uchar_t *ext_cigar, const uchar_t* tmp_ref_sequence, const uchar_t* tmp_read_sequence, ref_pos_t pos, uint32_t seq_len,
	uint32_t left_clipping, uint32_t right_clipping, uint32_t left_match, uint32_t right_match, int32_t del_size, 
	const scoring_t &scoring, double& affine_score)
{
	const char *decode = "ACGTNNXX";	
	
	if (left_match == 0 && right_match == 0) {
		std::fill(ext_cigar, ext_cigar + left_clipping, '$');
		std::fill(ext_cigar + left_clipping, ext_cigar + seq_len - right_clipping, '.');
		std::fill(ext_cigar + seq_len - right_clipping, ext_cigar + seq_len, '$');
		ext_cigar[seq_len] = 0;

		/*
		affine_score = (seq_len - left_clipping - right_clipping) * scoring.match +
			((left_clipping > 0) ? scoring.clipping : 0) +
			((right_clipping > 0) ? scoring.clipping : 0);
		*/

		affine_score = (seq_len - left_clipping - right_clipping) * scoring.match +
			((left_clipping > 0) ? scoring.gap_open + (left_clipping - 1) * scoring.gap_extend : 0) +
			((right_clipping > 0) ? scoring.gap_open + (right_clipping - 1) * scoring.gap_extend : 0);

	} else {
		
		std::fill(ext_cigar, ext_cigar + left_match, '.');
		
		int cigar_pos = left_match;
		int read_pos = left_match;
		pos += left_match;
		
		if (del_size > 0) {
			// deletion wrt reference - write corresponding reference symbols
			for (int i = 0; i < del_size; ++i, ++pos) {
				ext_cigar[cigar_pos++] = '^';
				ext_cigar[cigar_pos++] = decode[tmp_ref_sequence[pos]];
			}
		}
		else {
			// insertion wrt reference - write corresponding read symbols
			

			for (int i = 0; i < -del_size; ++i, ++read_pos) {
				ext_cigar[cigar_pos++] = '#';
				ext_cigar[cigar_pos++] = decode[tmp_read_sequence[read_pos + 1]];
			}
		}

		affine_score = (left_match + right_match) * scoring.match +
			scoring.gap_open + (abs(del_size) - 1) * scoring.gap_extend;

		std::fill(ext_cigar + cigar_pos, ext_cigar + cigar_pos + right_match, '.');
		cigar_pos += right_match;
		ext_cigar[cigar_pos] = 0;
	}
}

// ************************************************************************************
void CSoftClipping::clipIllegalPositions(mapping_desc_t& mapping, const scoring_t& scoring, const std::vector<seq_desc_t>& seq_desc, CMemoryPool<uchar_t>& mp_ext_cigar) {

	auto& ref_desc = seq_desc[mapping.ref_seq_id];
	int32_t ref_size = ref_desc.size + ref_desc.no_final_Ns + ref_desc.no_initial_Ns;
	int32_t mapping_last = mapping.ref_seq_pos + mapping.ref_length - 1;

	int toClipLeft = 0;
	int toClipRight = 0;

	if (mapping.ref_seq_pos <= 0) {
		toClipLeft = -mapping.ref_seq_pos + 1;
	}

	if (mapping_last > ref_size) {
		toClipRight = mapping_last - ref_size;
	}

	if (toClipLeft || toClipRight) {
		// fixme: perfect matching outside chromosome - this shouldn't happen
		if (mapping.method == MatchingMethod::Perfect || mapping.method == MatchingMethod::PerfectDistant) {
			mp_ext_cigar.Reserve(mapping.ext_cigar);
			std::fill(mapping.ext_cigar, mapping.ext_cigar + mapping.read_length, '.');
			mapping.ext_cigar[mapping.read_length] = 0;
		}
		
		mapping.left_clipping += toClipLeft;
		mapping.right_clipping += toClipRight;
		mapping.err_edit_distance -= toClipLeft + toClipRight;

		mapping.ref_length -= toClipLeft + toClipRight;
		mapping.ref_seq_pos += toClipLeft;

		mapping.score -= (toClipLeft + toClipRight) * (scoring.error + scoring.match);

		// only mismatches can be present as boundary mappings
		uchar_t * in = mapping.ext_cigar;
		while (toClipLeft) {
			*in++ = '$';
			--toClipLeft;
		}

		if (toClipRight) {
			in = mapping.ext_cigar + std::strlen(reinterpret_cast<char*>(mapping.ext_cigar)) - toClipRight;
			while (toClipRight) {
				*in++ = '$';
				--toClipRight;
			}
		}
	}
}

// ************************************************************************************
void CSoftClipping::clipLowQualityPositions(mapping_desc_t& mapping, const uchar_t* quality) {
	
	if (mapping.err_edit_distance == 0)
		return;

	uint32_t max_clipping = mapping.read_length - 20;

	string_reader_t quality_reader(quality, mapping.read_length, mapping.dir);

	// left side
	uchar_t* in = mapping.ext_cigar + mapping.left_clipping;
	uchar_t * out = in;

	while ((quality_reader[mapping.left_clipping] <= '#') && (((uint32_t) (mapping.left_clipping + mapping.right_clipping)) <= max_clipping)) {
		if (*in == '.') {
			// match - clip position, move in read and reference
			*out++ = '$';
			++in;
			++mapping.left_clipping;
			++mapping.ref_seq_pos;
			--mapping.ref_length;
		}
		else if (*in == '#') {
			// insertion - clip position, move in read, decrese edit distance
			*out++ = '$';
			in += 2;
			++mapping.left_clipping;
			--mapping.err_edit_distance;
		}
		else if (*in == '^') {
			// deletion - omit position, move in reference, decrese edit distance
			in += 2;
			++mapping.ref_seq_pos;
			--mapping.ref_length;
			--mapping.err_edit_distance;
		}
		else if (isalpha(*in)) {
			// mismatch - clip position, move in read and reference, decrese edit distance
			*out++ = '$';
			++in;
			++mapping.left_clipping;
			++mapping.ref_seq_pos;
			--mapping.ref_length;
			--mapping.err_edit_distance;
		}
	}

	if (in != out) {
		std::strcpy((char*)out, (char*)in); // out is on the left
	}

	// right side
	in = mapping.ext_cigar + std::strlen(reinterpret_cast<char*>(mapping.ext_cigar)) - mapping.right_clipping - 1;
	out = in;

	while ((quality_reader[mapping.read_length - mapping.right_clipping - 1] <= '#') && (mapping.left_clipping + mapping.right_clipping <= max_clipping)) {
		if (*in == '.') {
			// match - clip position, move in read and reference
			*out-- = '$';
			--in;
			++mapping.right_clipping;
			--mapping.ref_length;
		}
		else if (isalpha(*in)) {
			if (*(in - 1) == '#') {
				// insertion - clip position, move in read, decrese edit distance
				*out-- = '$';
				in -= 2;
				++mapping.right_clipping;
				--mapping.err_edit_distance;
			}
			else if (*(in - 1) == '^') {
				// deletion - omit position, move in reference, decrese edit distance
				in -= 2;
				--mapping.ref_length;
				--mapping.err_edit_distance;
			}
			else {

				// mismatch - clip position, move in read and reference, decrese edit distance
				*out-- = '$';
				--in;
				++mapping.right_clipping;
				--mapping.ref_length;
				--mapping.err_edit_distance;
			}
		}
	}

	if (in != out) {
		std::strcpy((char*)in, (char*)out); // in is on the left
	}
}

// ************************************************************************************
void  CSoftClipping::clipBoundaryPositions(mapping_desc_t& mapping, const scoring_t& scoring) 
{		
	if (mapping.err_edit_distance == 0)
		return;
	
	std::pair<int, double> mini{ -1, std::numeric_limits<double>::max() };
	std::pair<int, double> maxi{ -1, -std::numeric_limits<double>::max() };

	double cur = 0;

	// consider clipping penalty only once
	double leftPenalty = mapping.left_clipping > 0 ? 0 : scoring.clipping;
	double rightPenalty = mapping.right_clipping > 0 ? 0 : scoring.clipping;

	// left side
	uchar_t* in = mapping.ext_cigar + mapping.left_clipping;
	uchar_t * out = in;
	bool inIndel = false;

	int readPos = 0;
	int refPos = 0;

	for (int i = 0; i < std::max(mapping.read_length, mapping.ref_length) ; ++i) {
		if (*in == '.') {
			// match 
			cur += scoring.match;
			inIndel = false;
			++in;
			++refPos;
			++readPos;
		}
		else if (*in == '#') {
			// insertion wrt reference
			cur += inIndel ? scoring.gap_extend : scoring.gap_open;
			inIndel = true;
			in += 2;
			++readPos;
		}
		else if (*in == '^') {
			// deletion wrt reference
			cur += inIndel ? scoring.gap_extend : scoring.gap_open;
			inIndel = true;
			in += 2;
			++refPos;
		}
		else if (isalpha(*in)) {
			// mismatch
			cur += scoring.mismatch;
			inIndel = false;
			++in;
			++refPos;
			++readPos;
		}

		if (cur < mini.second) {
			mini.first = i;
			mini.second = cur;
		}

		if (cur > maxi.second) {
			maxi.first = i;
			maxi.second = cur;
		}
	}

	// left clipping
	if (mini.second + leftPenalty > 0) {

		in = mapping.ext_cigar + mapping.left_clipping;

		for (int i = 0; i < mini.first + 1; ++i) { // clip up to mini position
			if (*in == '.') {
				// match - clip position, move in read and reference
				*out++ = '$';
				++in;
				++mapping.left_clipping;
				++mapping.ref_seq_pos;
				--mapping.ref_length;
		
			}
			else if (*in == '#') {
				// insertion - clip position, move in read, decrese edit distance
				*out++ = '$';
				in += 2;
				++mapping.left_clipping;
				--mapping.err_edit_distance;
			
			}
			else if (*in == '^') {
				// deletion - omit position, move in reference, decrese edit distance
				in += 2;
				++mapping.ref_seq_pos;
				--mapping.ref_length;
				--mapping.err_edit_distance;
			
			}
			else if (isalpha(*in)) {
				// mismatch - clip position, move in read and reference, decrese edit distance
				*out++ = '$';
				++in;
				++mapping.left_clipping;
				++mapping.ref_seq_pos;
				--mapping.ref_length;
				--mapping.err_edit_distance;
			}
		}

		while (*in) {
			*out++ = *in++;
		}
		*out = 0;
	}

	// right clipping
	if (maxi.second + rightPenalty > cur) {

		in = mapping.ext_cigar + std::strlen(reinterpret_cast<char*>(mapping.ext_cigar)) - mapping.right_clipping - 1;
		out = in + 1;

	
		for (int i = 0; i < mapping.read_length - maxi.first; ++i) {
			if (*in == '.') {
				// match - clip position, move in read and reference
				*(--out) = '$';
				--in;
				++mapping.right_clipping;
				--mapping.ref_length;
			}
			else if (isalpha(*in)) {
				if (*(in - 1) == '#') {
					// insertion - clip position, move in read, decrese edit distance
					*(--out) = '$';
					in -= 2;
					++mapping.right_clipping;
					--mapping.err_edit_distance;
				}
				else if (*(in - 1) == '^') {
					// deletion - omit position, move in reference, decrese edit distance
					in -= 2;
					--mapping.ref_length;
					--mapping.err_edit_distance;
				}
				else {

					// mismatch - clip position, move in read and reference, decrese edit distance
					*(--out) = '$';
					--in;
					++mapping.right_clipping;
					--mapping.ref_length;
					--mapping.err_edit_distance;
				}
			}
		}

		++in;  // in is on the left
		while (*in) {
			*in++ = *out++;
		}
		*in = 0;
	}
}

// ************************************************************************************
void CSoftClipping::clipOverlappingPairs(mapping_pair_t& p, const scoring_t& scoring, CMemoryPool<uchar_t>& mp_ext_cigar) 
{
	if (p.first->ref_seq_id != p.second->ref_seq_id 
		|| p.first->left_clipping || p.first->right_clipping 
		|| p.second->left_clipping || p.second->right_clipping) {
		return;
	}
	
	mapping_desc_t* left = nullptr;
	mapping_desc_t* right = nullptr;

	if (p.first->ref_seq_pos < p.second->ref_seq_pos) {
		left = p.first;
		right = p.second;
	}
	else if (p.first->ref_seq_pos == p.second->ref_seq_pos) {
		if (p.first->dir == genome_t::direct) { // if positions are equal - treat direc as left
			left = p.first;
			right = p.second;
		}
		else {
			left = p.second;
			right = p.first;
		}
	}
	else {
		left = p.second;
		right = p.first;
	}

	// if reveresed read is on the left and overlaps with forward read
	if (left->dir == genome_t::rev_comp && left->ref_seq_pos + left->ref_length > right->ref_seq_pos) {
		for (int r = 0; r < 2; ++r) {
			auto& mapping = p[r];
			if (mapping.method == MatchingMethod::Perfect || mapping.method == MatchingMethod::PerfectDistant) {
				mp_ext_cigar.Reserve(mapping.ext_cigar);
				std::fill(mapping.ext_cigar, mapping.ext_cigar + mapping.read_length, '.');
				mapping.ext_cigar[mapping.read_length] = 0;
			}
		}
		
		int toClipLeft = right->ref_seq_pos - left->ref_seq_pos;
		int toClipRight = (right->ref_seq_pos + right->ref_length) - (left->ref_seq_pos + left->ref_length);

		clipLeft(*left, toClipLeft, scoring);
		clipRight(*right, toClipRight, scoring);
	}
}

// ************************************************************************************
void CSoftClipping::clipLeft(mapping_desc_t& mapping, int symbolsCount, const scoring_t& scoring) 
{
	uchar_t *in = mapping.ext_cigar + mapping.left_clipping;
	uchar_t * out = in;

	int clipped = 0;
	while (clipped < symbolsCount || (clipped == symbolsCount && *in=='#')) {
		if (*in == '.') {
			// match - clip position, move in read and reference
			*out++ = '$';
			++in;
			++mapping.left_clipping;
			++mapping.ref_seq_pos;
			--mapping.ref_length;
			++clipped;

		}
		else if (*in == '#') {
			// insertion - clip position, move in read, decrese edit distance
			*out++ = '$';
			in += 2;
			++mapping.left_clipping;
			--mapping.err_edit_distance;
		}
		else if (*in == '^') {
			// deletion - omit position, move in reference, decrese edit distance
			in += 2;
			++mapping.ref_seq_pos;
			--mapping.ref_length;
			--mapping.err_edit_distance;
			++clipped;

		}
		else if (isalpha(*in)) {
			// mismatch - clip position, move in read and reference, decrese edit distance
			*out++ = '$';
			++in;
			++mapping.left_clipping;
			++mapping.ref_seq_pos;
			--mapping.ref_length;
			--mapping.err_edit_distance;
			++clipped;
		}
	}

	if (out != in) {
		while (*in) {
			*out++ = *in++;
		}
		*out = 0;
	}
}

// ************************************************************************************
void CSoftClipping::clipRight(mapping_desc_t& mapping, int symbolsCount, const scoring_t& scoring) 
{
	uchar_t* in = mapping.ext_cigar + std::strlen(reinterpret_cast<char*>(mapping.ext_cigar)) - mapping.right_clipping - 1;
	uchar_t* out = in + 1;

	int clipped= 0;
	while (clipped < symbolsCount || (clipped == symbolsCount && *(in - 1) == '#')) {
		if (*in == '.') {
			// match - clip position, move in read and reference
			*(--out) = '$';
			--in;
			++mapping.right_clipping;
			--mapping.ref_length;
			++clipped;
		}
		else if (isalpha(*in)) {
			if (*(in - 1) == '#') {
				// insertion - clip position, move in read, decrese edit distance
				*(--out) = '$';
				in -= 2;
				++mapping.right_clipping;
				--mapping.err_edit_distance;
			}
			else if (*(in - 1) == '^') {
				// deletion - omit position, move in reference, decrese edit distance
				in -= 2;
				--mapping.ref_length;
				--mapping.err_edit_distance;
				++clipped;
			}
			else {

				// mismatch - clip position, move in read and reference, decrese edit distance
				*(--out) = '$';
				--in;
				++mapping.right_clipping;
				--mapping.ref_length;
				--mapping.err_edit_distance;
				++clipped;
			}
		}
	}

	++in;  // in is on the left
	if (in != out) {
		while (*in) {
			*in++ = *out++;
		}
		*in = 0;
	}
}

// EOF
