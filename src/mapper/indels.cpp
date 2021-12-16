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
#ifdef _MSC_VER
#pragma warning(disable : 6255)
#endif

#include "indels.h"

#define GET_CHAR_FROM_PATTERN(i, pattern_ptr) ((i) & 1 ? LO_NIBBLE(pattern_ptr[(i)>>1]) : HI_NIBBLE(pattern_ptr[(i)>>1]))
#define GET_CHAR_FROM_GENOME(j) ((left_end_index_in_gen + (j)) & 1 ? LO_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]) : HI_NIBBLE(ref_ptr[(left_end_index_in_gen+(j))>>1]))

// *******************************************************************************************
CIndelMatching::CIndelMatching(CParams *params)
{
//	gap_open = params->gap_open;
//	gap_extend = params->gap_extend;
	gap_ins_open = params->gap_ins_open;
	gap_ins_extend = params->gap_ins_extend;
	gap_del_open = params->gap_del_open;
	gap_del_extend = params->gap_del_extend;
	mismatch_score = params->mismatch_score;
	match_score = params->match_score;
	clipping_score = params->clipping_score;
	min_clipped_factor = params->min_clipped_factor;

	ref_ptr = nullptr;
	ref_size = 0;
	cur_offset = 0;
	query_len = 0;
	genome_prefetch = nullptr;

	rev_comp_code[(int32_t)sym_code_A] = sym_code_T;
	rev_comp_code[(int32_t)sym_code_C] = sym_code_G;
	rev_comp_code[(int32_t)sym_code_G] = sym_code_C;
	rev_comp_code[(int32_t)sym_code_T] = sym_code_A;
	rev_comp_code[4] = 4;
	rev_comp_code[5] = 5;
	rev_comp_code[6] = 6;
	rev_comp_code[7] = 7;

	fill_n(rev_comp_code + 8, 8, 15);

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

	fill_n(code2symbol, 16, 'N');

	code2symbol[(int32_t)sym_code_A] = 'A';
	code2symbol[(int32_t)sym_code_C] = 'C';
	code2symbol[(int32_t)sym_code_G] = 'G';
	code2symbol[(int32_t)sym_code_T] = 'T';
	code2symbol[(int32_t)sym_code_N_read] = 'N';
	code2symbol[(int32_t)sym_code_N_ref] = 'N';

	instruction_set_t instruction_set;
	int x = instrset_detect();
	if (x >= 0 && x <= 8)
		instruction_set = (instruction_set_t)x;
	else if (x < 0)
		instruction_set = instruction_set_t::none;
	else
		instruction_set = instruction_set_t::avx2;

	switch (instruction_set) {
	case instruction_set_t::none:
	case instruction_set_t::sse:
		throw new std::runtime_error("SSE2 extensions required!");
	case instruction_set_t::sse2:
	case instruction_set_t::sse3:
	case instruction_set_t::sse3s:
	case instruction_set_t::sse41:
	case instruction_set_t::sse42:
		ptr_CountMismatches = CountMismatches128<instruction_set_t::sse2>;
		ptr_PrefetchDecompressed = PrefetchDecompressed128<instruction_set_t::sse2>;
		break;
	case instruction_set_t::avx:
		ptr_CountMismatches = CountMismatches128<instruction_set_t::avx>;
		ptr_PrefetchDecompressed = PrefetchDecompressed128<instruction_set_t::avx>;
		break;
	case instruction_set_t::avx2:
		ptr_CountMismatches = CountMismatches256<instruction_set_t::avx2>;
		ptr_PrefetchDecompressed = PrefetchDecompressed256<instruction_set_t::avx2>;
		break;
	}
}

// *******************************************************************************************
CIndelMatching::~CIndelMatching()
{
}

// *******************************************************************************************
void CIndelMatching::SetReference(uchar_t* _ref_ptr, uint32_t _ref_size, uint32_t _cur_offset)
{
	ref_ptr = _ref_ptr;
	ref_size = _ref_size;
	cur_offset = _cur_offset;
}

// *******************************************************************************************
//
bool CIndelMatching::Match(ref_pos_t ref_pos, uint32_t max_indel_size, uint32_t fixed_segment_pos, uint32_t fixed_segment_size, uint32_t max_no_mismatches,
	candidate_mapping_t &candidate_mapping)
{
	candidate_mapping.no_mismatches = max_no_mismatches + 1;
//	candidate_mapping.type = mapping_type_t::none;
	uint32_t text_size = query_len + 2 * max_indel_size;

	uint32_t c_no_mismatches;
	uint32_t indel_pos;
	double min_acceptable_score = 0;
	double total_cost = min_acceptable_score;
	int min_acceptable_size_after_clipping = query_len * min_clipped_factor;

	// Left part data
	int left_indel_size = 0;
	uint32_t no_left_mismatches_indel = query_len;
	double best_left_cost = 0;
	uint32_t best_left_indel_pos = 0;
	bool left_indel_found = false;
	double best_left_clipping_cost = 0;
	uint32_t best_no_left_clipped = 0;
	uint32_t no_left_mismatches_clipped = 0;
	uint32_t best_no_left_mismatches_clipped = 0;

	// Right part data
	int right_indel_size = 0;
	uint32_t no_right_mismatches_indel = query_len;
	double best_right_cost = 0;
	uint32_t best_right_indel_pos = 0;
	bool right_indel_found = false;
	double best_right_clipping_cost = 0;
	uint32_t best_no_right_clipped = 0;
	uint32_t no_right_mismatches_clipped = 0;
	uint32_t best_no_right_mismatches_clipped = 0;

	v_left_candidates.clear();
	v_right_candidates.clear();

	prefetch_genome(ref_pos - max_indel_size, text_size);

	// **********************************
	// Verificatgion of fixed segment
	for (uint32_t i = 0; i < fixed_segment_size; ++i)
		if (query[fixed_segment_pos + i] != genome[max_indel_size + fixed_segment_pos + i])
			return false;

	// **********************************
	// Evaluate left part
	uint32_t no_left_mismatches_raw = 0;
	double clipping_cost = clipping_score;
	for (int i = (int)fixed_segment_pos - 1; i >= 0; --i)
	{
		if (query[i] != genome[max_indel_size + i])
		{
			++no_left_mismatches_raw;
			++no_left_mismatches_clipped;
			clipping_cost += mismatch_score;
		}
		else
		{
			clipping_cost += match_score;

			if (clipping_cost > best_left_clipping_cost)
			{
				best_left_clipping_cost = clipping_cost;
				best_no_left_clipped = i;
				best_no_left_mismatches_clipped = no_left_mismatches_clipped;
			}
		}
	}

	bcsw_query.Reset(min_side_size);

	if (fixed_segment_pos > (uint32_t) min_side_size)
	{
		for (int i = 0; i < min_side_size; ++i)
			bcsw_query.Insert(query[i], min_side_size - i - 1);

		find_candidates(bcsw_query, v_left_candidates, 0, max_indel_size + fixed_segment_pos);

		for (auto genome_left_pos : v_left_candidates)
		{
			tie(c_no_mismatches, indel_pos) = calc_mismatches_with_indel(genome_left_pos, max_indel_size + fixed_segment_pos + 1, 0, fixed_segment_pos + 1, true);
			if (c_no_mismatches > max_no_mismatches)
				continue;

			int c_indel = (int)max_indel_size - (int)genome_left_pos;

			if (abs(c_indel) > (int) max_indel_size)
				continue;

			double act_left_cost = (fixed_segment_pos - c_no_mismatches) * match_score + c_no_mismatches * mismatch_score;
			if (c_indel > 0)
				act_left_cost += gap_del_open + (c_indel - 1) * gap_del_extend;
			else if (c_indel < 0)
				act_left_cost += gap_ins_open + (-c_indel - 1) * gap_ins_extend;

			if (act_left_cost > best_left_cost)
			{
				best_left_cost = act_left_cost;
				no_left_mismatches_indel = c_no_mismatches;
				best_left_indel_pos = indel_pos;
				left_indel_size = c_indel;
				left_indel_found = c_indel != 0;
			}
		}
	}

	// **********************************
	// Evaluate right part
	uint32_t no_right_mismatches_raw = 0;
	clipping_cost = clipping_score;
	for (uint32_t i = fixed_segment_pos + fixed_segment_size; i < query_len; ++i)
	{
		if (query[i] != genome[max_indel_size + i])
		{
			++no_right_mismatches_raw;
			++no_right_mismatches_clipped;
			clipping_cost += mismatch_score;
		}
		else
		{
			clipping_cost += match_score;

			if (clipping_cost > best_right_clipping_cost)
			{
				best_right_clipping_cost = clipping_cost;
				best_no_right_mismatches_clipped = no_right_mismatches_clipped;
				best_no_right_clipped = query_len - i - 1;
			}
		}
	}

	if (fixed_segment_pos + fixed_segment_size + min_side_size < query_len)
	{
		bcsw_query.Clear();
		for (int i = 0; i < min_side_size; ++i)
			bcsw_query.Insert(query[(int64_t)query_len - min_side_size + i], min_side_size - i - 1);

		find_candidates(bcsw_query, v_right_candidates, max_indel_size + fixed_segment_pos + fixed_segment_size, text_size);

		for (auto right_pos : v_right_candidates)
		{
			tie(c_no_mismatches, indel_pos) = calc_mismatches_with_indel(fixed_segment_pos + fixed_segment_size - 1 + max_indel_size, right_pos + min_side_size,
				fixed_segment_pos + fixed_segment_size - 1, query_len, false);
			if (c_no_mismatches > max_no_mismatches)
				continue;

			int c_indel = ((int)right_pos + (int)min_side_size) - ((int)text_size - (int) max_indel_size);

			if (abs(c_indel) > (int) max_indel_size)
				continue;

			double act_right_cost = (query_len - fixed_segment_pos - fixed_segment_size - c_no_mismatches) * match_score + c_no_mismatches * mismatch_score;
			if (c_indel > 0)
				act_right_cost += gap_del_open + (c_indel - 1) * gap_del_extend;
			else if (c_indel < 0)
				act_right_cost += gap_ins_open + (-c_indel - 1) * gap_ins_extend;

			if (act_right_cost > best_right_cost)
			{
				best_right_cost = act_right_cost;
				no_right_mismatches_indel = c_no_mismatches;
				best_right_indel_pos = indel_pos;
				right_indel_size = c_indel;
				right_indel_found = c_indel != 0;
			}
		}
	}

	// **********************************
	// Case: indel - mismatches
	if (left_indel_found && no_left_mismatches_indel + no_right_mismatches_raw < max_no_mismatches)
	{
		double act_cost = best_left_cost + (query_len - fixed_segment_pos - no_right_mismatches_raw) * match_score + no_right_mismatches_raw * mismatch_score;
		if (act_cost > total_cost)
		{
			total_cost = act_cost;

			if(left_indel_size > 0)
				candidate_mapping.pos = ref_pos - (uint32_t) left_indel_size;
			else
				candidate_mapping.pos = ref_pos + (uint32_t)-left_indel_size;
			candidate_mapping.type = mapping_type_t::indel1;
			candidate_mapping.matching_len1 = best_left_indel_pos;
			candidate_mapping.indel_len1 = left_indel_size;
			candidate_mapping.no_mismatches = no_left_mismatches_indel + no_right_mismatches_raw;
			candidate_mapping.calc_penalty();
		}
	}
		
	// **********************************
	// Case: mismatches - indel
	if (right_indel_found && no_left_mismatches_raw + no_right_mismatches_indel < max_no_mismatches)
	{
		double act_cost = best_right_cost + (fixed_segment_pos + fixed_segment_size - no_left_mismatches_raw) * match_score + no_left_mismatches_raw * mismatch_score;
		if (act_cost > total_cost)
		{
			total_cost = act_cost;

			candidate_mapping.pos = ref_pos;
			candidate_mapping.type = mapping_type_t::indel1;
			candidate_mapping.matching_len1 = best_right_indel_pos;
			candidate_mapping.indel_len1 = right_indel_size;
			candidate_mapping.no_mismatches = no_left_mismatches_raw + no_right_mismatches_indel;
			candidate_mapping.calc_penalty();
		}
	}

	// **********************************
	// Case: indel - indel
	if (left_indel_found && right_indel_found && no_left_mismatches_indel + no_right_mismatches_indel < max_no_mismatches)
	{
		double act_cost = best_left_cost + best_right_cost + fixed_segment_size * match_score;
		if (act_cost > total_cost)
		{
			total_cost = act_cost;

			candidate_mapping.type = mapping_type_t::indel2;
			if (left_indel_size > 0)
				candidate_mapping.pos = ref_pos - (uint32_t)left_indel_size;
			else
				candidate_mapping.pos = ref_pos + (uint32_t)-left_indel_size;
			candidate_mapping.matching_len1 = best_left_indel_pos;
			candidate_mapping.indel_len1 = left_indel_size;

			if (left_indel_size > 0)
				candidate_mapping.matching_len2 = best_right_indel_pos - best_left_indel_pos;
			else
				candidate_mapping.matching_len2 = best_right_indel_pos - (best_left_indel_pos + (uint32_t) (-left_indel_size));

			candidate_mapping.indel_len2 = right_indel_size;
			candidate_mapping.no_mismatches = no_left_mismatches_indel + no_right_mismatches_indel;
			candidate_mapping.calc_penalty();
		}
	}

	// **********************************
	// Case: clipping - mismatches
	if (best_no_left_clipped > 0 && best_no_left_mismatches_clipped + no_right_mismatches_raw < rescale_max_no_mismatches(query_len - best_no_left_clipped, max_no_mismatches))
	{
		double act_cost = best_left_clipping_cost + (query_len - fixed_segment_pos - no_right_mismatches_raw) * match_score + no_right_mismatches_raw * mismatch_score;
		if (act_cost > total_cost && query_len - best_no_left_clipped >= min_acceptable_size_after_clipping)
		{
			total_cost = act_cost;

			candidate_mapping.pos = ref_pos + best_no_left_clipped;
			candidate_mapping.type = mapping_type_t::clipping_mismatches;
			candidate_mapping.clipping_len = best_no_left_clipped;
			candidate_mapping.matching_len1 = query_len - best_no_left_clipped;
			candidate_mapping.no_mismatches = best_no_left_mismatches_clipped + no_right_mismatches_raw;
			candidate_mapping.calc_penalty();
		}
	}

	// **********************************
	// Case: mismatches - clipping
	if (best_no_right_clipped > 0 && no_left_mismatches_raw + best_no_right_mismatches_clipped < rescale_max_no_mismatches(query_len - best_no_right_clipped, max_no_mismatches))
	{
		double act_cost = best_right_clipping_cost + (fixed_segment_pos + fixed_segment_size - no_left_mismatches_raw) * match_score + no_left_mismatches_raw * mismatch_score;
		if (act_cost > total_cost && query_len - best_no_right_clipped >= min_acceptable_size_after_clipping)
		{
			total_cost = act_cost;

			candidate_mapping.pos = ref_pos;
			candidate_mapping.type = mapping_type_t::mismatches_clipping;
			candidate_mapping.clipping_len = best_no_right_clipped;
			candidate_mapping.matching_len1 = query_len - best_no_right_clipped;
			candidate_mapping.no_mismatches = no_left_mismatches_raw + best_no_right_mismatches_clipped;
			candidate_mapping.calc_penalty();
		}
	}

	// **********************************
	// Case: clipping - indel
	if (best_no_left_clipped > 0 && right_indel_found && best_no_left_mismatches_clipped + no_right_mismatches_indel < rescale_max_no_mismatches(query_len - best_no_left_clipped, max_no_mismatches))
	{
		double act_cost = best_left_clipping_cost + best_right_cost + fixed_segment_size * match_score;
		if (act_cost > total_cost  && query_len - best_no_left_clipped >= min_acceptable_size_after_clipping)
		{
			total_cost = act_cost;

			candidate_mapping.pos = ref_pos + best_no_left_clipped;
			candidate_mapping.type = mapping_type_t::clipping_indel;
			candidate_mapping.indel_len1 = right_indel_size;
			candidate_mapping.clipping_len = best_no_left_clipped;
			candidate_mapping.matching_len1 = best_right_indel_pos - best_no_left_clipped;
			candidate_mapping.no_mismatches = best_no_left_mismatches_clipped + no_right_mismatches_indel;
			candidate_mapping.calc_penalty();
		}
	}	

	// **********************************
	// Case: indel - clipping
	if (best_no_right_clipped > 0 && left_indel_found && best_no_right_mismatches_clipped + no_left_mismatches_indel < rescale_max_no_mismatches(query_len - best_no_right_clipped, max_no_mismatches))
	{
		double act_cost = best_right_clipping_cost + best_left_cost + fixed_segment_size * match_score;
		if (act_cost > total_cost && query_len - best_no_right_clipped >= min_acceptable_size_after_clipping)
		{
			total_cost = act_cost;

			if (left_indel_size > 0)
				candidate_mapping.pos = ref_pos - (uint32_t)left_indel_size;
			else
				candidate_mapping.pos = ref_pos + (uint32_t)-left_indel_size;

			candidate_mapping.type = mapping_type_t::indel_clipping;
			candidate_mapping.indel_len1 = left_indel_size;
			candidate_mapping.clipping_len = best_no_right_clipped;
			candidate_mapping.matching_len1 = best_left_indel_pos;
			candidate_mapping.no_mismatches = best_no_right_mismatches_clipped + no_left_mismatches_indel;
			candidate_mapping.calc_penalty();
		}
	}

	// **********************************
	// Case: clipping - clipping
	if (best_no_left_clipped > 0 && best_no_right_clipped > 0 && 
		best_no_left_mismatches_clipped + best_no_right_mismatches_clipped < rescale_max_no_mismatches(query_len - best_no_left_clipped - best_no_right_clipped, max_no_mismatches))
	{
		double act_cost = best_left_clipping_cost + best_right_clipping_cost + fixed_segment_size * match_score;
		if (act_cost > total_cost && query_len - best_no_left_clipped - best_no_right_clipped >= min_acceptable_size_after_clipping)
		{
			total_cost = act_cost;

			candidate_mapping.pos = ref_pos + best_no_left_clipped;
			candidate_mapping.type = mapping_type_t::clipping_clipping;
			candidate_mapping.clipping_len = best_no_left_clipped;
			candidate_mapping.matching_len1 = query_len - best_no_left_clipped - best_no_right_clipped;
			candidate_mapping.no_mismatches = best_no_left_mismatches_clipped + best_no_right_mismatches_clipped;
			candidate_mapping.calc_penalty();
		}
	}


	return candidate_mapping.type != mapping_type_t::none;
}

// ************************************************************************************
void CIndelMatching::prefetch_genome(uint32_t begin, uint32_t size)
{
	genome.resize((size_t)size + 2ll);

	(*ptr_PrefetchDecompressed)(genome.data(), ref_ptr + (begin >> 1), size, begin & 1);
}

// ************************************************************************************
bool CIndelMatching::find_candidates(CBitCharsSingleWord &bcsw_query, vector<uint32_t>& v_cand, uint32_t begin, uint32_t end)
{
	uint32_t min_no_matches = min_side_size - max_no_side_errors;

	v_cand.clear();
	
	bcsw_ref.Reset(min_side_size);

	for (int i = 0; i < min_side_size - 1; ++i)
		bcsw_ref.Insert(genome[(int64_t) begin + i], min_side_size - i - 1 - 1);

	for (int i = begin + min_side_size - 1; i < (int) end; ++i)
	{
		bcsw_ref.Insert(genome[i]);

		if (bcsw_ref.CountMatches(bcsw_query) >= min_no_matches)
			v_cand.emplace_back(i - min_side_size + 1);
	}

	return !v_cand.empty();
}

// ************************************************************************************
pair<uint32_t, uint32_t> CIndelMatching::calc_mismatches_with_indel(uint32_t genome_left_pos, uint32_t genome_right_pos, uint32_t query_left_pos, uint32_t query_right_pos, bool left_margin)
{
	uint32_t text_size = genome_right_pos - genome_left_pos;
	uint32_t query_size = query_right_pos - query_left_pos;
	uint32_t stats_size = min(text_size, query_size);
	
	if (stats_size <= (uint32_t) min_side_size)
		return make_pair(query_len + 1, 0);
	
	int start_pos = left_margin ? min_side_size : 0;
	int end_pos = stats_size - (left_margin ? 0: min_side_size);

	//int del_size = (int)text_size - (int)query_size;

	v_left_mismatches.resize(stats_size + 1ull);
	v_right_mismatches.resize(stats_size + 1ull);

	uchar_t *p = &genome[genome_left_pos];
	uchar_t* q = &query[query_left_pos];
	auto dest = v_left_mismatches.begin();
	*dest++ = 0;

	uint32_t prev_val = 0;
	for (uint32_t i = 1; i <= stats_size; ++i)
	{
		*dest = prev_val + (*p++ != *q++);
		prev_val = *dest++;
	}

	p = &genome[(int64_t)genome_left_pos + text_size - 1];
	q = &query[query_left_pos + (int64_t)query_size - 1];
	dest = v_right_mismatches.begin() + stats_size;
	*dest-- = 0;

	prev_val = 0;
	for (uint32_t i = 1; i <= stats_size; ++i)
	{
		*dest = prev_val + (*p-- != *q--);
		prev_val = *dest--;
	}

	int best_no_mismatches = query_len + 1;
	int best_indel_pos = 0;

	for (int i = start_pos; i < end_pos; ++i)
	{
		int tot_mis = v_left_mismatches[i] + v_right_mismatches[i];
		if (tot_mis < best_no_mismatches)
		{
			best_indel_pos = i;
			best_no_mismatches = tot_mis;
		}
	}

	return make_pair(best_no_mismatches, query_left_pos + best_indel_pos);
}

// ************************************************************************************
bool CIndelMatching::GetExtCigar(ref_pos_t ref_pos, int32_t del_size, uchar_t* ext_cigar, uint32_t& no_mismatches)
{
	int text_size = (int)query_len + del_size;
	genome_prefetch = (uchar_t*)alloca(text_size + 10ll);

	// Genome prefetch
	uchar_t* __restrict ptr = genome_prefetch;
	uchar_t* genome_ptr = ref_ptr + (ref_pos >> 1);
	uint32_t text_size_div2 = (text_size + 3) / 2;

	if (ref_pos & 1)
		*ptr++ = *genome_ptr++ & 0x0f;

	for (uint32_t i = 0; i < text_size_div2; ++i)
	{
		*ptr++ = *genome_ptr >> 4;
		*ptr++ = *genome_ptr++ & 0x0f;
	}

	int stats_size = min(text_size, (int)query_len);

	int *mismatches_left = (int*) alloca((stats_size + 1ll) * sizeof(int));
	int *mismatches_right = (int*) alloca((stats_size + 1ll) * sizeof(int));

	mismatches_left[0] = 0;
	for (int i = 1; i <= stats_size; ++i)
		mismatches_left[i] = mismatches_left[i - 1ll] + (genome_prefetch[i - 1ll] != query[i - 1ll]);

	mismatches_right[stats_size] = 0;
	for (int i = 1; i <= stats_size; ++i)
		mismatches_right[(int64_t) stats_size - i] = mismatches_right[(int64_t) stats_size - i + 1ll] + (genome_prefetch[(int64_t) text_size - i] != query[(int64_t) query_len - i]);

	int best_no_mismatches = query_len + 1;
	int left_part_size = 0;
	int right_part_size = 0;

	if (del_size > 0)
	{
		for (int i = min_side_size; i < stats_size - min_side_size; ++i)
		{
			int tot_mis = mismatches_left[i] + mismatches_right[i];
			if (tot_mis < best_no_mismatches) 
			{
				best_no_mismatches = tot_mis;
				left_part_size = i;
				right_part_size = query_len - i;
			}
		}
	}
	else
	{
		for (int i = min_side_size; i < stats_size - min_side_size; ++i)
		{
			int tot_mis = mismatches_left[i] + mismatches_right[i];
			if (tot_mis < best_no_mismatches)
			{
				best_no_mismatches = tot_mis;
				left_part_size = i;
				right_part_size = text_size - i;
			}
		}
	}

	no_mismatches = best_no_mismatches;

	uchar_t* p_ec = ext_cigar;

	for (int i = 0; i < left_part_size; ++i)
	{
		if (genome_prefetch[i] == query[i])
			*p_ec++ = '.';
		else
			*p_ec++ = code2symbol[genome_prefetch[i]];
	}

	int i_ref, i_query;

	if (del_size > 0)
	{
		for (int i = 0; i < del_size; ++i)
		{
			*p_ec++ = '^';
			*p_ec++ = code2symbol[genome_prefetch[left_part_size + i]];
		}

		i_ref = left_part_size + del_size;
		i_query = left_part_size;
	}
	else
	{
		for (int i = 0; i < -del_size; ++i)
		{
			*p_ec++ = '#';
			*p_ec++ = code2symbol[query[(int64_t) left_part_size + i]];
		}

		i_ref = left_part_size;
		i_query = left_part_size - del_size;
	}

	for (int i = 0; i < right_part_size; ++i)
	{
		if (genome_prefetch[i_ref + i] == query[(int64_t) i_query + i])
			* p_ec++ = '.';
		else
			*p_ec++ = code2symbol[genome_prefetch[i_ref + i]];
	}

	*p_ec = 0;

	return true;
}

// ************************************************************************************
void CIndelMatching::append_matching_part(uint32_t mp_len, uchar_t* &ec, uint32_t &i_genome, uint32_t &i_query, 
	uint16_t& err_edit_distance, uint16_t& num_events, double& score, CParams* params)
{
	int no_mis = 0;

	for (uint32_t i = 0; i < mp_len; ++i)
	{
		if (genome[i_genome] == query[i_query])
			*ec++ = '.';
		else
		{
			*ec++ = code2symbol[genome[i_genome]];
			++no_mis;
		}

		++i_genome;
		++i_query;
	}

	score += (mp_len - no_mis) * params->match_score + no_mis * params->mismatch_score;
	num_events += no_mis;
	err_edit_distance += no_mis;
}

// ************************************************************************************
void CIndelMatching::append_indel(int indel_len, uchar_t* &ec, uint32_t& i_genome, uint32_t& i_query, 
	uint16_t &err_edit_distance, uint16_t & num_events, double& score, CParams* params)
{
	if (indel_len > 0)
	{
		for (uint32_t i = 0; i < (uint32_t)indel_len; ++i)
		{
			*ec++ = '^';
			*ec++ = code2symbol[genome[i_genome++]];
		}

		score += params->gap_del_open + (indel_len - 1) * params->gap_del_extend;
		num_events += 1;
		err_edit_distance += indel_len;
	}
	else
	{
		for (uint32_t i = 0; i < (uint32_t)-indel_len; ++i)
		{
			*ec++ = '#';
			*ec++ = code2symbol[query[i_query++]];
		}

		score += params->gap_ins_open + (indel_len - 1) * params->gap_ins_extend;
		num_events += 1;
		err_edit_distance += -indel_len;
	}
}

// ************************************************************************************
void CIndelMatching::append_clipping(int clipping_len, uchar_t*& ec, uint32_t& i_genome, uint32_t& i_query,
	uint16_t& err_edit_distance, uint16_t& num_events, double& score, CParams* params)
{
	for (uint32_t i = 0; i < (uint32_t) clipping_len; ++i)
		*ec++ = '$';

	score += params->clipping_score;
	num_events += 1;

	i_query += clipping_len;
}

// ************************************************************************************
bool CIndelMatching::get_ext_cigar_indel1(mapping_desc_t& mapping_desc, CParams* params)
{
	uint32_t i_genome = 0;
	uint32_t i_query = 0;
	uchar_t* ec = mapping_desc.ext_cigar;

	append_matching_part(mapping_desc.mapping.matching_len1, ec, i_genome, i_query, 
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);
	append_indel(mapping_desc.mapping.indel_len1, ec, i_genome, i_query, 
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);
	append_matching_part(mapping_desc.read_length - i_query, ec, i_genome, i_query, 
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

	*ec = 0;

	if(mapping_desc.mapping.indel_len1 > 0)
		mapping_desc.method = MatchingMethod::Enumeration::MappingDelete;
	else
		mapping_desc.method = MatchingMethod::Enumeration::MappingInsert;

	mapping_desc.ref_length = (uint16_t) ((int)mapping_desc.read_length + (int)mapping_desc.mapping.indel_len1);
	mapping_desc.left_clipping = 0;
	mapping_desc.right_clipping = 0;

	return true;
}

// ************************************************************************************
bool CIndelMatching::get_ext_cigar_indel2(mapping_desc_t& mapping_desc, CParams* params)
{
	uint32_t i_genome = 0;
	uint32_t i_query = 0;
	uchar_t* ec = mapping_desc.ext_cigar;

	append_matching_part(mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);
	append_indel(mapping_desc.mapping.indel_len1, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);
	append_matching_part(mapping_desc.mapping.matching_len2, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);
	append_indel(mapping_desc.mapping.indel_len2, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);
	append_matching_part(mapping_desc.read_length - i_query, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

	*ec = 0;

	if (mapping_desc.mapping.indel_len1 > 0)
		if(mapping_desc.mapping.indel_len2 > 0)
			mapping_desc.method = MatchingMethod::Enumeration::MappingDeleteDelete;
		else
			mapping_desc.method = MatchingMethod::Enumeration::MappingDeleteInsert;
	else
		if (mapping_desc.mapping.indel_len2 > 0)
			mapping_desc.method = MatchingMethod::Enumeration::MappingInsertDelete;
		else
			mapping_desc.method = MatchingMethod::Enumeration::MappingInsertInsert;

	mapping_desc.ref_length = (uint16_t)((int)mapping_desc.read_length + (int)mapping_desc.mapping.indel_len1 + (int)mapping_desc.mapping.indel_len2);
	mapping_desc.left_clipping = 0;
	mapping_desc.right_clipping = 0;

	return true;
}

// ************************************************************************************
bool CIndelMatching::get_ext_cigar_clipping_mismatches(mapping_desc_t& mapping_desc, CParams* params)
{
	uint32_t i_genome = 0;
	uint32_t i_query = 0;
	uchar_t* ec = mapping_desc.ext_cigar;

	if (mapping_desc.mapping.type == mapping_type_t::clipping_mismatches)
	{
		mapping_desc.method = MatchingMethod::Enumeration::MappingClippingMismatches;

		append_clipping(mapping_desc.read_length - mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		mapping_desc.left_clipping = mapping_desc.read_length - mapping_desc.mapping.matching_len1;
	}

	append_matching_part(mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

	if (mapping_desc.mapping.type == mapping_type_t::mismatches_clipping)
	{
		mapping_desc.method = MatchingMethod::Enumeration::MappingMismatchesClipping;

		append_clipping(mapping_desc.read_length - mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		mapping_desc.right_clipping = mapping_desc.read_length - mapping_desc.mapping.matching_len1;
	}

	*ec = 0;

	mapping_desc.ref_length = (uint16_t)((int)mapping_desc.mapping.matching_len1);

	return true;
}

// ************************************************************************************
bool CIndelMatching::get_ext_cigar_clipping_indel(mapping_desc_t& mapping_desc, CParams* params)
{
	// !!! nieskonczone
	uint32_t i_genome = 0;
	uint32_t i_query = 0;
	uchar_t* ec = mapping_desc.ext_cigar;

	if (mapping_desc.mapping.type == mapping_type_t::clipping_indel)
	{
		if(mapping_desc.mapping.indel_len1 > 0)
			mapping_desc.method = MatchingMethod::Enumeration::MappingClippingDelete;
		else
			mapping_desc.method = MatchingMethod::Enumeration::MappingClippingInsert;

		append_clipping(mapping_desc.mapping.clipping_len, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		append_matching_part(mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		append_indel(mapping_desc.mapping.indel_len1, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		append_matching_part(query_len - i_query, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		mapping_desc.left_clipping = mapping_desc.mapping.clipping_len;
	}
	else
	{
		if (mapping_desc.mapping.indel_len1 > 0)
			mapping_desc.method = MatchingMethod::Enumeration::MappingDeleteClipping;
		else
			mapping_desc.method = MatchingMethod::Enumeration::MappingInsertClipping;

		append_matching_part(mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		append_indel(mapping_desc.mapping.indel_len1, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		append_matching_part(query_len - i_query - mapping_desc.mapping.clipping_len, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		append_clipping(mapping_desc.mapping.clipping_len, ec, i_genome, i_query,
			mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

		mapping_desc.right_clipping = mapping_desc.mapping.clipping_len;
	}

	*ec = 0;

	mapping_desc.ref_length = (uint16_t)((int)mapping_desc.read_length - (int)mapping_desc.mapping.clipping_len + (int) mapping_desc.mapping.indel_len1);

	return true;
}

// ************************************************************************************
bool CIndelMatching::get_ext_cigar_clipping_clipping(mapping_desc_t& mapping_desc, CParams* params)
{
	uint32_t i_genome = 0;
	uint32_t i_query = 0;
	uchar_t* ec = mapping_desc.ext_cigar;

	mapping_desc.method = MatchingMethod::Enumeration::MappingClippingClipping;

	mapping_desc.left_clipping = mapping_desc.mapping.clipping_len;

	append_clipping(mapping_desc.left_clipping, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

	append_matching_part(mapping_desc.mapping.matching_len1, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

	append_clipping(mapping_desc.read_length - mapping_desc.mapping.matching_len1 - mapping_desc.left_clipping, ec, i_genome, i_query,
		mapping_desc.err_edit_distance, mapping_desc.num_events, mapping_desc.score, params);

	mapping_desc.right_clipping = mapping_desc.read_length - mapping_desc.mapping.matching_len1 - mapping_desc.left_clipping;

	*ec = 0;

	mapping_desc.ref_length = (uint16_t)((int)mapping_desc.mapping.matching_len1);

	return true;
}

// ************************************************************************************
bool CIndelMatching::GetExtCigarNew(mapping_desc_t& mapping_desc, CParams *params)
{
	mapping_desc.score = 0.0;
	mapping_desc.err_edit_distance = 0;
	mapping_desc.num_events = 0;
	mapping_desc.score_type = ScoringType::Affine;

	if (mapping_desc.mapping.type == mapping_type_t::indel1)
		return get_ext_cigar_indel1(mapping_desc, params);
	if (mapping_desc.mapping.type == mapping_type_t::indel2)
		return get_ext_cigar_indel2(mapping_desc, params);
	if(mapping_desc.mapping.type == mapping_type_t::clipping_mismatches || mapping_desc.mapping.type == mapping_type_t::mismatches_clipping)
		return get_ext_cigar_clipping_mismatches(mapping_desc, params);
	if(mapping_desc.mapping.type == mapping_type_t::clipping_indel || mapping_desc.mapping.type == mapping_type_t::indel_clipping)
		return get_ext_cigar_clipping_indel(mapping_desc, params);
	if(mapping_desc.mapping.type == mapping_type_t::clipping_clipping)
		return get_ext_cigar_clipping_clipping(mapping_desc, params);

	return false;
}

// ************************************************************************************
bool CIndelMatching::Preprocess(uchar_t* seq, uint32_t seq_len, genome_t orientation)
{
	query_len = seq_len;

	query.resize(query_len);

	if (orientation == genome_t::direct)
		for (uint32_t i = 0; i < seq_len; ++i)
			query[i] = (uchar_t)GET_CHAR_FROM_PATTERN(i, seq);
	else
		for (uint32_t i = 0; i < seq_len; ++i)
			query[i] = rev_comp_code[GET_CHAR_FROM_PATTERN(seq_len - i - 1, seq)];

	return true;
}

// ************************************************************************************
bool CIndelMatching::PreprocessRawSeq(uchar_t* seq, uint32_t seq_len, genome_t orientation)
{
	query_len = seq_len;

	query.resize(query_len);

	if (orientation == genome_t::direct)
		for (uint32_t i = 0; i < seq_len; ++i)
			query[i] = raw_code[seq[i]];
	else
		for (uint32_t i = 0; i < seq_len; ++i)
			query[i] = raw_rev_comp_code[seq[seq_len - i - 1]];

	return true;
}

// ************************************************************************************
uint32_t CIndelMatching::indel_to_mismatch_score(int32_t size)
{
	double affine_score = 0;
	
	if(size > 0)
		affine_score = gap_del_open + ((double) size - 1) * gap_del_extend;
	else
		affine_score = gap_ins_open + ((double)-size - 1) * gap_ins_extend;

	return (uint32_t)(affine_score / mismatch_score);
}

// EOF
