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

#include "sam.h"

#pragma warning (disable: 6255 6263)

#ifdef ENABLE_VCF_VARIANTS
const uint32_t min_indel_len = 3u;
const uint32_t max_indel_len = 50u;

// ************************************************************************************
template <>
bool CSamGenerator::find_pairs<MatchingMethod::Enumeration::VariantInsertLong>(
	std::vector<mapping_desc_t> mapping_desc[2],
	std::vector<mapping_desc_t*> hits[2],
	uchar_t *id[2],
	uchar_t *sequence[2],
	uint32_t sequence_len[2],
	uchar_t *quality[2],
	std::vector<mapping_pair_t>& mapping_pairs,
	std::pair<int, int>& editDistanceRange)
{
	string paired_ref_seq_name;
	int32_t paired_ref_seq_pos;
	int32_t paired_ref_seq_id;

	Evaluator <mapping_pair_t> pairEvaluator(*insertSizeModel);

	// use model if enough samples were collected
	ref_pos_t distanceThreshold = (ref_pos_t) (insertSizeModel->mean() + insertSizeModel->getVariantDev());

	vector<uint32_t> v_pat_pos;

	for (uint32_t r = 0; r < 2; ++r) {
		for (auto h : hits[r]) {
/*			if (h->err_edit_distance < editDistanceRange.first || h->err_edit_distance > editDistanceRange.second || h->method == MatchingMethod::LevMyers) {
				continue; // ignore hits if edit distance is outside specified range
			}*/
			
			mapping_desc_t& this_hit = *h;
			genome_t paired_dir = get_direction_using_orientation(this_hit.mapping.direction, mapping_orientation);

			// Prepare sequences for comparison
			// Note: 0th positions in tmp_read_sequence and tmp_ref_sequences are guards (just to make the DP calculation easier)
			ref_pos_t start_pos =
				(paired_dir == genome_t::direct) ?
				(this_hit.mapping.pos - distanceThreshold + sequence_len[!r]) :
				(this_hit.mapping.pos + 1);

			auto var_range_ins = variant_db->FindRangeIns(start_pos, start_pos + distanceThreshold);

			if (paired_dir == genome_t::direct) {
				copy_direct(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
			}
			else {
				convert_to_rev_comp(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
			}
//			ref_copy_direct(tmp_ref_sequence, reference->GetData(), start_pos, distanceThreshold);

			const int min_pattern_len = 10;
			for (auto p = var_range_ins.first; p != var_range_ins.second; ++p)
			{
				uint8_t *q = variant_db->GetInsDesc(p->ins_desc_pos);
				int len = p->ins_len;
				uchar_t *pattern;
				int pattern_len;

				if (len < min_indel_len)
					continue;		// !!! tymczasowo nie sprawdzam krotkich insercji

				if (len > max_indel_len)
					continue;		// !!! tymczasowo nie sprawdzam d³ugich insercji

				if (len < min_pattern_len)
				{
					pattern = (uchar_t*) alloca(min_pattern_len);
					uint32_t ref_pos = p->pos;
					auto ref_ptr = reference->GetData();
					for (int i = 1; i < min_pattern_len - len + 1; ++i)
						pattern[i + len - 1] = GET_SYMBOL(ref_ptr, ref_pos + i);

					pattern_len = min_pattern_len;
				}
				else
				{
					pattern = (uchar_t*)alloca(len);
					pattern_len = len;
				}

				for (int i = 0; i < len; ++i)
					pattern[i] = q[i];

				if (!variantMatch->find_pattern(tmp_read_sequence+1, sequence_len[!r], pattern, pattern_len, v_pat_pos))
					continue;

				// !!! Probuje znalezc matche na pozostalym kawalku
				for (auto cand_pos : v_pat_pos)
				{
					bool is_match = true;
					uint32_t no_mis = 0;
					uint32_t max_no_mis = (uint32_t)(sequence_len[!r] * max_mate_edit_distance_frac);
					auto ref_ptr = reference->GetData();

					for (uint32_t i = 0; i <= cand_pos && no_mis <= max_no_mis; ++i)
						if (GET_SYMBOL(ref_ptr, p->pos - i) != tmp_read_sequence[cand_pos - i + 1])
							++no_mis;

					for (uint32_t i = cand_pos + p->ins_len; i < sequence_len[!r] && no_mis <= max_no_mis; ++i)
						if (GET_SYMBOL(ref_ptr, p->pos + 1u + (i - (cand_pos + p->ins_len))) != tmp_read_sequence[i + 1])
							++no_mis;

					if(no_mis <= max_no_mis)
						cout << "ins found: " + to_string(p->pos) + "   of len: " + to_string(p->ins_len) + "    no_mis: " + to_string(no_mis) + "   read_pos: " + to_string(cand_pos) + "\n";

					if (no_mis <= max_no_mis)
					{
						uint32_t paired_raw_pos = p->pos - cand_pos;
						if (ref_seq_desc->Translate(paired_raw_pos, paired_ref_seq_name, paired_ref_seq_pos, paired_ref_seq_id, id[!r]))
						{
							uchar_t* ext_cigar;
							mp_ext_cigar->Reserve(ext_cigar);

							// Ustawianie Ext CIGAR
							uint32_t ext_cigar_pos = 0;

							for (uint32_t i = 0; i <= cand_pos; ++i)
								if (GET_SYMBOL(ref_ptr, p->pos - i) == tmp_read_sequence[cand_pos - i + 1])
									ext_cigar[cand_pos - ext_cigar_pos++] = '.';
								else
									ext_cigar[cand_pos - ext_cigar_pos++] =	code2dna[tmp_read_sequence[cand_pos - i + 1]];

							for (uint32_t i = 1; i < p->ins_len; ++i)
							{
								ext_cigar[ext_cigar_pos++] = '#';
								ext_cigar[ext_cigar_pos++] = code2dna[pattern[i]];
							}

							for (uint32_t i = cand_pos + p->ins_len; i < sequence_len[!r]; ++i)
								if (GET_SYMBOL(ref_ptr, p->pos + 1u + (uint32_t)(i - (cand_pos + p->ins_len))) == tmp_read_sequence[i + 1])
									ext_cigar[ext_cigar_pos++] = '.';
								else
									ext_cigar[ext_cigar_pos++] = code2dna[tmp_read_sequence[i + 1]];
							ext_cigar[ext_cigar_pos] = 0;

							if (mapping_desc[!r].size() == mapping_desc[!r].capacity())
								cerr << "Will restruct 1\n";

							candidate_mapping_t paired_mapping = candidate_mapping_t::construct<mapping_type_t::lev>(paired_raw_pos, paired_dir, 0);

							mapping_desc[!r].emplace_back(paired_mapping, paired_ref_seq_id, paired_ref_seq_pos, sequence_len[!r], ext_cigar,
								params->gap_del_open, params->gap_del_extend, params->gap_ins_open, params->gap_ins_extend, params->mismatch_score);
							auto& mate_hit = mapping_desc[!r].back();
							mate_hit.method = MatchingMethod::Enumeration::VariantInsertLong;
							mate_hit.ref_length = sequence_len[!r] + p->ins_len - 1;
							mate_hit.err_edit_distance = no_mis + p->ins_len - 1u;
							mate_hit.score = 
								((double) sequence_len[!r] - (double) no_mis - ((double) p->ins_len - 1.0)) * params->match_score + 
								no_mis * params->mismatch_score +
								((double) p->ins_len - 1 - 1) * params->gap_ins_extend + params->gap_ins_open;
							mate_hit.score_type = ScoringType::Affine;

							if (this_hit.method == MatchingMethod::Enumeration::Unmatched) {
								fill_cigar_with_lev(this_hit, id[r], sequence[r], sequence_len[r], quality[r], false, true);
							}

							softClipping->clipIllegalPositions(this_hit, insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);
							softClipping->clipIllegalPositions(mate_hit, insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);
							if (params->enable_boundary_clipping) {
								softClipping->clipBoundaryPositions(this_hit, insertSizeModel->getScoring());
								softClipping->clipBoundaryPositions(mate_hit, insertSizeModel->getScoring());
							}

							mapping_pair_t candidate(
								r == 0 ? &this_hit : &mate_hit,
								r == 0 ? &mate_hit : &this_hit);

							candidate.score = pairEvaluator(candidate);

							mapping_pairs.push_back(candidate);

#ifdef COLLECT_STATS
							running_stats->AddTotals((uint32_t) STAT_MAPPED_PE_HISTO_VAR_INS_LONG + min(STAT_MAX_LEN_INDELS, p->ins_len), (int64_t) 1);
#endif
						}
					}
				}
			}
		}
	}

	if (mapping_pairs.size()) {
//		std::sort_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);
		std::sort(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);

		// fill cigars in the mapping pair elements
		for (auto& mp : mapping_pairs) {
			for (int r = 0; r < 2; ++r) {
				if (mp[r].method == MatchingMethod::Enumeration::VariantInsertLong) {
					hits[r].push_back(&mp[r]);
				}
			}
		}

//		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][0].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
//		++stat_mapped_pe_histo_clipped[min((size_t)mapping_pairs[0][1].err_edit_distance / 5, stat_mapped_pe_histo_clipped.size() - 1)];
	}

	return mapping_pairs.size() > 0;
}


// ************************************************************************************
template <>
bool CSamGenerator::find_pairs<MatchingMethod::Enumeration::VariantDeleteLong>(
	std::vector<mapping_desc_t> mapping_desc[2],
	std::vector<mapping_desc_t*> hits[2],
	uchar_t *id[2],
	uchar_t *sequence[2],
	uint32_t sequence_len[2],
	uchar_t *quality[2],
	std::vector<mapping_pair_t>& mapping_pairs,
	std::pair<int, int>& editDistanceRange)
{
	string paired_ref_seq_name;
	int32_t paired_ref_seq_pos;
	int32_t paired_ref_seq_id;

	Evaluator <mapping_pair_t> pairEvaluator(*insertSizeModel);

	// use model if enough samples were collected
	ref_pos_t distanceThreshold = (ref_pos_t) (insertSizeModel->mean() + insertSizeModel->getVariantDev());

	vector<uint32_t> v_pat_pos;

	for (uint32_t r = 0; r < 2; ++r) {
		for (auto h : hits[r]) {
			/*			if (h->err_edit_distance < editDistanceRange.first || h->err_edit_distance > editDistanceRange.second || h->method == MatchingMethod::LevMyers) {
							continue; // ignore hits if edit distance is outside specified range
						}*/

			mapping_desc_t& this_hit = *h;
			genome_t paired_dir = get_direction_using_orientation(this_hit.mapping.direction, mapping_orientation);

			// Prepare sequences for comparison
			// Note: 0th positions in tmp_read_sequence and tmp_ref_sequences are guards (just to make the DP calculation easier)
			ref_pos_t start_pos =
				(paired_dir == genome_t::direct) ?
				(this_hit.mapping.pos - distanceThreshold + sequence_len[!r]) :
				(this_hit.mapping.pos + 1);

			auto var_range_del = variant_db->FindRangeDel(start_pos, start_pos + distanceThreshold);

			if (paired_dir == genome_t::direct) {
				copy_direct(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
			}
			else {
				convert_to_rev_comp(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
			}
			//			ref_copy_direct(tmp_ref_sequence, reference->GetData(), start_pos, distanceThreshold);

			const int min_pattern_len = 10;
			for (auto p = var_range_del.first; p != var_range_del.second; ++p)
			{
				int len = p->del_len;
				uchar_t *pattern;
				int pattern_len;

				if (len < min_indel_len)
					continue;		// !!! tymczasowo nie sprawdzam krotkich delecji
				if (len > max_indel_len)
					continue;		// !!! tymczasowo nie sprawdzam d³ugich insercji

				uint32_t ref_pos = p->pos;
				auto ref_ptr = reference->GetData();
				pattern_len = min_pattern_len;
				pattern = (uchar_t*)alloca(pattern_len);
				pattern[0] = GET_SYMBOL(ref_ptr, ref_pos);
				for (int i = 1; i < pattern_len; ++i)
					pattern[i] = GET_SYMBOL(ref_ptr, ref_pos + p->del_len - 1u + i);

				if (!variantMatch->find_pattern(tmp_read_sequence + 1, sequence_len[!r], pattern, pattern_len, v_pat_pos))
					continue;

				// !!! Probuje znalezc matche na pozostalym kawalku
				for (auto cand_pos : v_pat_pos)
				{
					bool is_match = true;
					uint32_t no_mis = 0;
					uint32_t max_no_mis = (uint32_t)(sequence_len[!r] * max_mate_edit_distance_frac);
					auto ref_ptr = reference->GetData();

					for (uint32_t i = 0; i <= cand_pos && no_mis <= max_no_mis; ++i)
						if (GET_SYMBOL(ref_ptr, p->pos - i) != tmp_read_sequence[cand_pos - i + 1u])
							++no_mis;

					for (uint32_t i = cand_pos + 1; i < sequence_len[!r] && no_mis <= max_no_mis; ++i)
						if (GET_SYMBOL(ref_ptr, p->pos + i - (cand_pos + 1u) + p->del_len) != tmp_read_sequence[i + 1])
							++no_mis;

					if (no_mis <= max_no_mis)
						cout << "del found: " + to_string(p->pos) + "   of len: " + to_string(p->del_len) + "    no_mis: " + to_string(no_mis) + "   read_pos: " + to_string(cand_pos) + "\n";

					if (no_mis <= max_no_mis)
					{
						uint32_t paired_raw_pos = p->pos - cand_pos;
						if (ref_seq_desc->Translate(paired_raw_pos, paired_ref_seq_name, paired_ref_seq_pos, paired_ref_seq_id, id[!r]))
						{
							uchar_t* ext_cigar;
							mp_ext_cigar->Reserve(ext_cigar);

							// Ustawianie Ext CIGAR
							uint32_t ext_cigar_pos = 0;

							for (uint32_t i = 0; i <= cand_pos; ++i)
								if (GET_SYMBOL(ref_ptr, p->pos - i) == tmp_read_sequence[cand_pos - i + 1u])
									ext_cigar[cand_pos - ext_cigar_pos++] = '.';
								else
									ext_cigar[cand_pos - ext_cigar_pos++] = code2dna[tmp_read_sequence[cand_pos - i + 1]];

							for (uint32_t i = 1; i < p->del_len; ++i)
							{
								ext_cigar[ext_cigar_pos++] = '^';
								ext_cigar[ext_cigar_pos++] = code2dna[GET_SYMBOL(ref_ptr, p->pos+i)];
							}

							for (uint32_t i = cand_pos + 1; i < sequence_len[!r] && no_mis <= max_no_mis; ++i)
								if (GET_SYMBOL(ref_ptr, p->pos + i - (cand_pos + 1u) + p->del_len) == tmp_read_sequence[i + 1])
									ext_cigar[ext_cigar_pos++] = '.';
								else
									ext_cigar[ext_cigar_pos++] = code2dna[tmp_read_sequence[i + 1]];
							ext_cigar[ext_cigar_pos] = 0;

							if (mapping_desc[!r].size() == mapping_desc[!r].capacity())
								cerr << "Will restruct 2\n";

							candidate_mapping_t paired_mapping = candidate_mapping_t::construct<mapping_type_t::lev>(paired_raw_pos, paired_dir, 0);

							mapping_desc[!r].emplace_back(paired_mapping, paired_ref_seq_id, paired_ref_seq_pos, sequence_len[!r], ext_cigar,
								params->gap_del_open, params->gap_del_extend, params->gap_ins_open, params->gap_ins_extend, params->mismatch_score);
							auto& mate_hit = mapping_desc[!r].back();
							mate_hit.method = MatchingMethod::Enumeration::VariantDeleteLong;
							mate_hit.ref_length = sequence_len[!r] - p->del_len + 1;
							mate_hit.err_edit_distance = no_mis + p->del_len - 1u;
							mate_hit.score =
								((double) sequence_len[!r] - no_mis) * params->match_score +
								no_mis * params->mismatch_score +
								((double) p->del_len - 1u - 1u) * params->gap_del_extend + params->gap_del_open;
							mate_hit.score_type = ScoringType::Affine;

							if (this_hit.method == MatchingMethod::Enumeration::Unmatched) {
								fill_cigar_with_lev(this_hit, id[r], sequence[r], sequence_len[r], quality[r], false, true);
							}

							softClipping->clipIllegalPositions(this_hit, insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);
							softClipping->clipIllegalPositions(mate_hit, insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);
							if (params->enable_boundary_clipping) {
								softClipping->clipBoundaryPositions(this_hit, insertSizeModel->getScoring());
								softClipping->clipBoundaryPositions(mate_hit, insertSizeModel->getScoring());
							}

							mapping_pair_t candidate(
								r == 0 ? &this_hit : &mate_hit,
								r == 0 ? &mate_hit : &this_hit);

							candidate.score = pairEvaluator(candidate);

							mapping_pairs.push_back(candidate);

#ifdef COLLECT_STATS
							running_stats->AddTotals((uint32_t)STAT_MAPPED_PE_HISTO_VAR_DEL_LONG + min(STAT_MAX_LEN_INDELS, p->del_len), (int64_t)1);
#endif
						}
					}
				}
			}
		}
	}

	if (mapping_pairs.size()) {
		//		std::sort_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);
		std::sort(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);

		// fill cigars in the mapping pair elements
		for (auto& mp : mapping_pairs) {
			for (int r = 0; r < 2; ++r) {
				if (mp[r].method == MatchingMethod::Enumeration::VariantDeleteLong) {
					hits[r].push_back(&mp[r]);
				}
			}
		}

		//		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][0].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
		//		++stat_mapped_pe_histo_clipped[min((size_t)mapping_pairs[0][1].err_edit_distance / 5, stat_mapped_pe_histo_clipped.size() - 1)];
	}

	return mapping_pairs.size() > 0;;
}

#endif

// EOF

