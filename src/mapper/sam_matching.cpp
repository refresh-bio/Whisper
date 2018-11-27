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


#include "sam.h"

//#define LEV_MYERS_ASSERTIONS

// ************************************************************************************
bool pushMappingHeap(std::vector<mapping_pair_t>& mapping_pairs, const mapping_pair_t& candidate, int heapSize) 
{
	bool stillBuilding = (int) mapping_pairs.size() < heapSize;
	bool replaceMinimum = ((int) mapping_pairs.size() == heapSize) && (candidate.score > mapping_pairs.front().score);
	
	if (stillBuilding || replaceMinimum) {

		if (stillBuilding) {
			mapping_pairs.push_back(candidate);
		}

		if (replaceMinimum) {
			pop_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);
			mapping_pairs.back() = candidate;
		}

		// min-heap - smallest score at the top of the heap
		push_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);
		return true;
	}

	return false;
}

// ************************************************************************************
void CSamGenerator::generate_unique_hits(std::vector<mapping_desc_t>& mapping_desc, std::vector<mapping_desc_t*>& hits) 
{
	hits.resize(mapping_desc.size());
	std::transform(mapping_desc.begin(), mapping_desc.end(), hits.begin(), [](mapping_desc_t& m)->mapping_desc_t* {
		return &m;
	});
	
	// Remove redundant mappings
	sort(hits.begin(), hits.end(), [](const mapping_desc_t* x, const mapping_desc_t* y) {
		if (x->ref_seq_id != y->ref_seq_id)
			return x->ref_seq_id < y->ref_seq_id;
		return x->ref_seq_pos < y->ref_seq_pos;
	});
	auto it = unique(hits.begin(), hits.end(), [](const mapping_desc_t *x, const mapping_desc_t *y) {
		return x->ref_seq_id == y->ref_seq_id && x->ref_seq_pos == y->ref_seq_pos;
	});

	hits.erase(it, hits.end());

	if (hits.size() > 1) {
		// iterate over all hits
		for (auto a = hits.begin(); a < std::prev(hits.end()); ++a) {
			auto best = a;
			auto prev_b = a;
			auto b = std::next(prev_b);
			// find groups of hits less than 
			for (b = std::next(a); b != hits.end(); ++b, ++prev_b) {
				auto ref_hit = (params->hit_merging_wrt_first) ? a : prev_b;

				if (
					((*b)->ref_seq_id != (*ref_hit)->ref_seq_id) ||
					((*b)->ref_seq_pos - (*ref_hit)->ref_seq_pos > (int) params->hit_merging_threshold)) {
					break;
				}
				
				if ((*b)->err_edit_distance < (*best)->err_edit_distance) {
					best = b;
				}
			}

			// replace current with best
			if (best != a) {
				*a = *best;
			}
			hits.erase(std::next(a), b); // does nothing if only one hit
		}
	}
}

// ************************************************************************************
bool CSamGenerator::find_duplicate_pairs(const std::vector<mapping_pair_t>& mapping_pairs) 
{
	auto tmp = mapping_pairs;
	std::sort(tmp.begin(), tmp.end(), mapping_pair_t::compareByPositions);

	static int cnt = 0;

	for (auto a = tmp.begin(); a < std::prev(tmp.end()); ++a) {
		for (auto b = std::next(a); b != tmp.end(); ++b) {

			if (b->first->ref_seq_pos - a->first->ref_seq_pos > 0) {
				break;
			}

			if (abs(b->second->ref_seq_pos - a->second->ref_seq_pos) == 0) {
				cerr << "Duplicate " << cnt++ << endl;
				//return true;
			}
		}
	}

	return false;
}

// ************************************************************************************
template <>
bool CSamGenerator::find_pairs<MatchingMethod::FromMapping>(
	std::vector<mapping_desc_t> mapping_desc[2],
	std::vector<mapping_desc_t*> hits[2],
	uchar_t *id[2],
	uchar_t *sequence[2],
	uint32_t sequence_len[2],
	uchar_t *quality[2],
	std::vector<mapping_pair_t>& mapping_pairs,
	std::pair<int, int>& editDistanceRange)
{
	if (hits[0].empty() || hits[1].empty()) {
		return false;
	}

	mapping_desc_t* bestHits[2];
	bestHits[0] = *std::min_element(hits[0].begin(), hits[0].end(), mapping_desc_t::compareByScoresDescending_ptr);
	bestHits[1] = *std::min_element(hits[1].begin(), hits[1].end(), mapping_desc_t::compareByScoresDescending_ptr);

	bool isSensitive[2];
	isSensitive[0] = (bestHits[0]->err_edit_distance > params->max_no_errors) ? true : false;
	isSensitive[1] = (bestHits[1]->err_edit_distance > params->max_no_errors) ? true : false;

	int bestId = (bestHits[0]->err_edit_distance < bestHits[1]->err_edit_distance) ? 0 : 1;

	if (isSensitive[bestId]) {
		return false;
	}


/*	if (isSensitive[0] || isSensitive[1]) {
		return false;
	}
*/
	Evaluator <mapping_pair_t> evaluator(*insertSizeModel);
	
	// Find pairs of mappings within range
	int lessHitsId = hits[0].size() > hits[1].size();

	std::pair<int32_t, int32_t> distanceRange{
		insertSizeModel->mean() - insertSizeModel->getHighConfidenceDev(), insertSizeModel->mean() + insertSizeModel->getHighConfidenceDev() };

	for (auto p : hits[lessHitsId]) {
		// create dummy hits just to establish search ranges 
		auto lo = *p;
		lo.ref_seq_pos -= distanceRange.second;
		auto hi = *p;
		hi.ref_seq_pos += distanceRange.second;
		
		auto begin = std::lower_bound(hits[!lessHitsId].begin(), hits[!lessHitsId].end(), &lo, mapping_desc_t::compareByPositions_ptr);
		auto end = std::upper_bound(hits[!lessHitsId].begin(), hits[!lessHitsId].end(), &hi, mapping_desc_t::compareByPositions_ptr);
		
		std::vector<mapping_desc_t*>& mate_desc = hits[!lessHitsId];
	//	for (auto q: hits[!lessHitsId]) {
		for (auto it = begin; it != mate_desc.end() && it != end; ++it) {
			auto& q = *it;
			int dist = mapping_pair_t::calculateInsertSize(*p, *q); // returns MAX_INT for different chromosomes

			if (dist >= distanceRange.first && dist <= distanceRange.second) {
				mapping_pair_t candidate(
					lessHitsId == 0 ? p : q,
					lessHitsId == 0 ? q : p);

				candidate.score = evaluator(candidate);

				pushMappingHeap(mapping_pairs, candidate, 2000);
			}
		}	
	}
		
	if (mapping_pairs.size()) {
		// sort mapping pairs decreasingly by score
		std::sort_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);

		if (mapping_pairs[0].first->err_edit_distance == bestHits[0]->err_edit_distance 
			&& mapping_pairs[0].second->err_edit_distance == bestHits[1]->err_edit_distance) {

			++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][0].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
			++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][1].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
		}
		else {
			mapping_pairs.resize(0);
		}
	}

	return mapping_pairs.size() > 0;
}

// ************************************************************************************
template <>
bool CSamGenerator::find_pairs<MatchingMethod::FromMappingDistant>(
	std::vector<mapping_desc_t> mapping_desc[2], 
	std::vector<mapping_desc_t*> hits[2],
	uchar_t *id[2], 
	uchar_t *sequence[2], 
	uint32_t sequence_len[2], 
	uchar_t *quality[2], 
	std::vector<mapping_pair_t>& mapping_pairs,
	std::pair<int, int>& editDistanceRange)
{
	if (hits[0].empty() || hits[1].empty()) {
		return false;
	}

	Evaluator <mapping_pair_t> evaluator(*insertSizeModel);

	// Find pairs of mappings within range
	int lessHitsId = hits[0].size() > hits[1].size();

	std::pair<int32_t, int32_t> distanceRanges[2]{
		{ insertSizeModel->mean() - insertSizeModel->getSaturationDev(), insertSizeModel->mean() - insertSizeModel->getMyersDev() },
		{ insertSizeModel->mean() + insertSizeModel->getMyersDev(), insertSizeModel->mean() + insertSizeModel->getSaturationDev() } };

	for (auto p : hits[lessHitsId]) {
		// create dummy hits just to establish search ranges 
		auto lo = *p;
		lo.ref_seq_pos -= distanceRanges[1].second;
		auto hi = *p;
		hi.ref_seq_pos += distanceRanges[1].second;

		std::vector<mapping_desc_t*>& mate_desc = hits[!lessHitsId];

		auto begin = std::lower_bound(mate_desc.begin(), mate_desc.end(), &lo, mapping_desc_t::compareByPositions_ptr);
		auto end = std::upper_bound(mate_desc.begin(), mate_desc.end(), &hi, mapping_desc_t::compareByPositions_ptr);

		for (auto it = begin; it != mate_desc.end() && it != end; ++it) {
			auto& q = *it;

			int dist = mapping_pair_t::calculateInsertSize(*p, *q); // returns MAX_INT for different chromosomes

			if ((dist >= distanceRanges[0].first && dist <= distanceRanges[0].second) || (dist >= distanceRanges[1].first && dist <= distanceRanges[1].second)) {
				mapping_pair_t candidate(
					lessHitsId == 0 ? p : q,
					lessHitsId == 0 ? q : p);

				candidate.score = evaluator(candidate);

				pushMappingHeap(mapping_pairs, candidate, 2000);
			}
		}
	}

	// fixme: take into account only single distant mapping
	mapping_desc_t* bestHits[2];
	bestHits[0] = *std::min_element(hits[0].begin(), hits[0].end(), mapping_desc_t::compareByScoresDescending_ptr);
	bestHits[1] = *std::min_element(hits[1].begin(), hits[1].end(), mapping_desc_t::compareByScoresDescending_ptr);

	mapping_pair_t distantMapping(bestHits[0], bestHits[1]);
	distantMapping.score = evaluator(distantMapping);

	// if outside saturation ranges
	if (distantMapping.calculateInsertSize() < distanceRanges[0].first || distantMapping.calculateInsertSize() > distanceRanges[1].second) {
		pushMappingHeap(mapping_pairs, distantMapping, 2000);
	}

	if (mapping_pairs.size()) {
		// sort mapping pairs decreasingly by score
		std::sort_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);
		
		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][0].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][1].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
	}

	return mapping_pairs.size() > 0;;
	
}

// ************************************************************************************
template <>
bool CSamGenerator::find_pairs<MatchingMethod::LevMyers>(
	std::vector<mapping_desc_t> mapping_desc[2],
	std::vector<mapping_desc_t*> hits[2],
	uchar_t *id[2],
	uchar_t *sequence[2],
	uint32_t sequence_len[2],
	uchar_t *quality[2],
	std::vector<mapping_pair_t>& mapping_pairs,
	std::pair<int, int>& editDistanceRange)
{
//	remove_repetitive(mapping_desc[0]);
//	remove_repetitive(mapping_desc[1]);

	LevMyers* levMyers[2];
	levMyers[0] = (sequence_len[0] < 256)
		? (sequence_len[0] < 128) ? levMyers128 : levMyers256
		: levMyers64;

	levMyers[1] = (sequence_len[1] < 256)
		? (sequence_len[1] < 128) ? levMyers128 : levMyers256
		: levMyers64;

	uint32_t paired_raw_pos;
	string paired_ref_seq_name;
	int32_t paired_ref_seq_pos;
	int32_t paired_ref_seq_id;

	Evaluator <mapping_pair_t> pairEvaluator(*insertSizeModel);

	// use model if enough samples were collected
	double distanceThreshold = insertSizeModel->mean() + insertSizeModel->getMyersDev();

	std::vector<mapping_desc_t> new_hits[2];
	
	double bestPairScores[2] = { -std::numeric_limits<double>::max() , -std::numeric_limits<double>::max() };

	for (uint32_t r = 0; r < 2; ++r) {
		// this is the highest mate score that can be found using LevMyers
		auto bestMate = std::min_element(hits[!r].begin(), hits[!r].end(), mapping_desc_t::compareByScoresDescending_ptr);

		double optimisticMateScore = (bestMate == hits[!r].end())
			? insertSizeModel->getScoring().getOptimisticScore(params->max_read_len, 0)
			: insertSizeModel->getScoring().getOptimisticScore((*bestMate)->read_length, (*bestMate)->err_edit_distance);
			
		for (auto h : hits[r]) {
			if (h->err_edit_distance < editDistanceRange.first || h->err_edit_distance > editDistanceRange.second) {
				continue; // ignore hits if edit distance is outside specified range
			}
			
			double optimisticThisScore = insertSizeModel->getScoring().getOptimisticScore(h->read_length, h->err_edit_distance);
			double allowedPairPenalty = -bestPairScores[1] + (optimisticThisScore + optimisticMateScore);
			
			// if there is no chance for pair to overcome second best
			if (allowedPairPenalty < insertSizeModel->getMinPenalty()) {
	//			continue;
			}

			// default TLEN range
	//		ref_pos_t minTlen = 0;
	//		ref_pos_t maxTlen = distanceThreshold;

			// narrow TLEN range according to the maximum pair penalty
			if (allowedPairPenalty < insertSizeModel->getMaxPenalty()) {
	//			minTlen = std::max(0, insertSizeModel->calculateInvPenalty(allowedPairPenalty));
	//			maxTlen = 2 * insertSizeModel->mean() - minTlen; // symmetric range
			}

			mapping_desc_t& this_hit = *h;
			genome_t paired_dir = get_direction_using_orientation(this_hit.dir, mapping_orientation);
			uint32_t max_mate_edit_distance = (uint32_t)(sequence_len[!r] * max_mate_edit_distance_frac) * 2;

			// Prepare sequences for comparison
			// Note: 0th positions in tmp_read_sequence and tmp_ref_sequences are guards (just to make the DP calculation easier)
			ref_pos_t start_pos = ((paired_dir == genome_t::direct) ? 
				(this_hit.raw_pos + sequence_len[!r] - distanceThreshold) :
				(this_hit.raw_pos + 1 - max_mate_edit_distance)) ;
			
			uint32_t edit_distance;
			ref_pos_t relative_end_pos;

			levMyers[!r]->preprocessRawSeq(sequence[!r], sequence_len[!r], paired_dir);
			bool ok = levMyers[!r]->dynamicProgramming(
				start_pos, 
				distanceThreshold + 2 * max_mate_edit_distance,
				max_mate_edit_distance, relative_end_pos, edit_distance);
		
#ifdef LEV_MYERS_ASSERTIONS			
			uint32_t pos2, edit_distance2;
			levMyers64->preprocessRawSeq(sequence[!r], sequence_len[!r], paired_dir);
			levMyers64->dynamicProgramming(start_pos, distanceThreshold, max_mate_edit_distance, pos2, edit_distance2);
			ASSERT((pos == pos2), " CSamGenerator::store_mapped_pair_reads() DP error (mappings = 0)!"); // fixme: why 1000?
#endif

			if (ok) {
				if (paired_dir == genome_t::direct) {
					copy_direct(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
				}
				else {
					convert_to_rev_comp(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
				}
				ref_copy_direct(tmp_ref_sequence, reference->GetData(), start_pos, distanceThreshold + 2 * max_mate_edit_distance);

				uchar_t* ext_cigar;
				double affine_score;
				uint32_t num_events;
				mp_ext_cigar->Reserve(ext_cigar);
				ref_pos_t relative_begin_pos = 
					levMyers[!r]->getExtCigar(ext_cigar, tmp_ref_sequence, tmp_read_sequence, quality[!r], 
						relative_end_pos, edit_distance, sequence_len[!r], paired_dir, insertSizeModel->getScoring(), affine_score, num_events);

#ifdef LEV_MYERS_ASSERTIONS		
				string q1((char*)ext_cigar);
				uchar_t ext_cigar2[1001];
				levMyers64->getExtCigar(ext_cigar2, tmp_ref_sequence, tmp_read_sequence, pos, edit_distance, sequence_len[!r]);
				string q2((char*)ext_cigar2);
				ASSERT(q1 == q2, "CSamGenerator::store_mapped_pair_reads() CIGAR extraction error (mappings = 0)");
#endif				
				paired_raw_pos = start_pos + relative_begin_pos;

				if (ref_seq_desc->Translate(paired_raw_pos, paired_ref_seq_name, paired_ref_seq_pos, paired_ref_seq_id, id[!r])) {
					//						if(((genome_t) paired_dir) == genome_t::rev_comp)
					//							paired_ref_seq_pos -= sequence_len[!r] - 1;
					// Mapping must be to the same reference sequence
					if (num_events <= max_mate_edit_distance / 2 && this_hit.ref_seq_id == paired_ref_seq_id) {
						
						if (this_hit.method == MatchingMethod::Unmatched) {
							fill_cigar_with_lev(this_hit, id[r], sequence[r], sequence_len[r], quality[r], false, true);
						}

						if (this_hit.num_events <= max_mate_edit_distance / 2) {

							mapping_desc[!r].emplace_back(paired_raw_pos, paired_ref_seq_id, paired_ref_seq_pos, paired_dir, sequence_len[!r], edit_distance, ext_cigar);
							auto& mate_hit = mapping_desc[!r].back();

							mate_hit.method = MatchingMethod::LevMyers;
							mate_hit.ref_length = relative_end_pos - relative_begin_pos;

							// make score affine
							mate_hit.score = affine_score;
							mate_hit.score_type = ScoringType::Affine;

							mapping_pair_t candidate(
								r == 0 ? &this_hit : &mate_hit,
								r == 0 ? &mate_hit : &this_hit);

							candidate.score = pairEvaluator(candidate);
							pushMappingHeap(mapping_pairs, candidate, 2000000);

							if (candidate.score > bestPairScores[0]) {
								bestPairScores[1] = bestPairScores[0];
								bestPairScores[0] = candidate.score;
							}
							else if (candidate.score < bestPairScores[0] && candidate.score > bestPairScores[1]) {
								bestPairScores[1] = candidate.score;
							}
						}
						else {
							mp_ext_cigar->Free(ext_cigar);
						}
					}
					else {
						mp_ext_cigar->Free(ext_cigar);
					}
				}
				else
				{
					cerr << "Myers: start_pos: " << start_pos << "   relative_end_pos: " << relative_end_pos << "   relative_begin_pos: " << relative_begin_pos << endl;
					cerr << ext_cigar << endl;
					mp_ext_cigar->Free(ext_cigar);
				}
			}
		}
	}

	if (mapping_pairs.size()) {

		// sort mapping pairs decreasingly by score
		std::sort_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);

		// fill cigars in the mapping pair elements
		for (auto& mp : mapping_pairs) {
			for (int r = 0; r < 2; ++r) {
				if (mp[r].method == MatchingMethod::LevMyers) {
					// add new hits in sorted list
					auto it = upper_bound(hits[r].begin(), hits[r].end(), &mp[r], mapping_desc_t::compareByPositions_ptr);
					hits[r].insert(it, &mp[r]);				
				} 
			}
		}
		
		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][0].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][1].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
	}
	
	return mapping_pairs.size() > 0;;
}

// ************************************************************************************
template <>
bool CSamGenerator::find_pairs<MatchingMethod::Clipping>(
	std::vector<mapping_desc_t> mapping_desc[2],
	std::vector<mapping_desc_t*> hits[2],
	uchar_t *id[2],
	uchar_t *sequence[2],
	uint32_t sequence_len[2],
	uchar_t *quality[2],
	std::vector<mapping_pair_t>& mapping_pairs,
	std::pair<int, int>& editDistanceRange)
{

	uint32_t paired_raw_pos;
	string paired_ref_seq_name;
	int32_t paired_ref_seq_pos;
	int32_t paired_ref_seq_id;

	Evaluator <mapping_pair_t> pairEvaluator(*insertSizeModel);

	// use model if enough samples were collected
	double distanceThreshold = insertSizeModel->mean() + insertSizeModel->getMyersDev();

	for (uint32_t r = 0; r < 2; ++r) {
		for (auto h : hits[r]) {
			if (h->err_edit_distance < editDistanceRange.first || h->err_edit_distance > editDistanceRange.second || h->method == MatchingMethod::LevMyers) {
				continue; // ignore hits if edit distance is outside specified range
			}
			
			
			mapping_desc_t& this_hit = *h;
			genome_t paired_dir = get_direction_using_orientation(this_hit.dir, mapping_orientation);
			
			// Prepare sequences for comparison
			// Note: 0th positions in tmp_read_sequence and tmp_ref_sequences are guards (just to make the DP calculation easier)
			ref_pos_t start_pos = 
				(paired_dir == genome_t::direct) ? 
				(this_hit.raw_pos - distanceThreshold + sequence_len[!r]) :
				(this_hit.raw_pos + 1);
			ref_pos_t relative_begin_pos;

			uint32_t left_clipping;
			uint32_t right_clipping;
			uint32_t left_match;
			uint32_t right_match;
			int32_t del_size;

			softClipping->preprocessRawSeq(sequence[!r], sequence_len[!r], paired_dir);
			if (softClipping->match(start_pos, distanceThreshold, HashTable::hash_len, relative_begin_pos, left_clipping, right_clipping,
				left_match, right_match, del_size)) {
				paired_raw_pos = start_pos + relative_begin_pos;
				if (ref_seq_desc->Translate(paired_raw_pos, paired_ref_seq_name, paired_ref_seq_pos, paired_ref_seq_id, id[!r])) {
					if (this_hit.ref_seq_id == paired_ref_seq_id) {
						if (paired_dir == genome_t::direct) {
							copy_direct(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
						}
						else {
							convert_to_rev_comp(tmp_read_sequence + 1, sequence[!r], sequence_len[!r]);
						}
						ref_copy_direct(tmp_ref_sequence, reference->GetData(), start_pos, distanceThreshold);
					
						uchar_t* ext_cigar;
						double affine_score;
						mp_ext_cigar->Reserve(ext_cigar);
						softClipping->getExtCigar(ext_cigar, tmp_ref_sequence, tmp_read_sequence, relative_begin_pos, sequence_len[!r], 
							left_clipping, right_clipping, left_match, right_match, del_size,
							insertSizeModel->getScoring(), affine_score);

						mapping_desc[!r].emplace_back(paired_raw_pos, paired_ref_seq_id, paired_ref_seq_pos, paired_dir, sequence_len[!r], 0, ext_cigar);
						auto& mate_hit = mapping_desc[!r].back();
						mate_hit.method = MatchingMethod::Clipping;
						mate_hit.ref_length = sequence_len[!r] - (left_clipping + right_clipping) + del_size;
						mate_hit.err_edit_distance = abs(del_size);
						mate_hit.left_clipping = left_clipping;
						mate_hit.right_clipping = right_clipping;
						mate_hit.ref_seq_pos += left_clipping;
						mate_hit.score = affine_score;
						mate_hit.score_type = ScoringType::Affine;

						if (this_hit.method == MatchingMethod::Unmatched) {
							fill_cigar_with_lev(this_hit, id[r], sequence[r], sequence_len[r], quality[r], false, true);
						}
						
						mapping_pair_t candidate(
							r == 0 ? &this_hit : &mate_hit,
							r == 0 ? &mate_hit : &this_hit);

						candidate.score = pairEvaluator(candidate);

						pushMappingHeap(mapping_pairs, candidate, 2000000);
					}
				} 
			}
		}
	}

	if (mapping_pairs.size()) {
		std::sort_heap(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresDescending);
	
		// fill cigars in the mapping pair elements
		for (auto& mp : mapping_pairs) {
			for (int r = 0; r < 2; ++r) {
				if (mp[r].method == MatchingMethod::Clipping) {
					// add new hits in sorted list
					auto it = upper_bound(hits[r].begin(), hits[r].end(), &mp[r], mapping_desc_t::compareByPositions_ptr);
					hits[r].insert(it, &mp[r]);
				} 
			}
		}

		++stat_mapped_pe_histo[min((size_t)mapping_pairs[0][0].err_edit_distance, stat_mapped_pe_histo.size() - 1)];
		++stat_mapped_pe_histo_clipped[min((size_t)mapping_pairs[0][1].err_edit_distance / 5, stat_mapped_pe_histo_clipped.size() - 1)];
	}

	return mapping_pairs.size() > 0;;
}

// ************************************************************************************
bool CSamGenerator::find_single(std::vector<mapping_desc_t*> hits[2], uchar_t *id[2], uchar_t *sequence[2], uint32_t sequence_len[2], uchar_t* quality[2]) {

	for (int i = 0; i < 2; ++i) {
		if (hits[i].size()) {
			auto& cur_hits = hits[i];
//			double affine_score;
			cur_hits.resize(1);
				
			// store only best alignment
			++stat_mapped_pe_histo[min((size_t)cur_hits[0]->err_edit_distance, stat_mapped_pe_histo.size() - 1)];

			fill_cigar_with_lev(*cur_hits[0], id[i], sequence[i], sequence_len[i], quality[i], false, false);
			softClipping->clipIllegalPositions(*cur_hits[0], insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);
		}
	}

	return true;
}

// ************************************************************************************
void CSamGenerator::fill_cigar_with_lev(mapping_desc_t& mapping, uchar_t *id, uchar_t *sequence, uint32_t sequence_len, uchar_t* quality, bool distant, bool set_affine_score)
{

	// perform LevMyers only when it makes sense (no perfect matching)
	if (mapping.err_edit_distance == 0) {
		mapping.ref_length = sequence_len;
		mapping.method = distant ? MatchingMethod::PerfectDistant : MatchingMethod::Perfect;
		mapping.num_events = 0;
		mapping.score = sequence_len * insertSizeModel->getScoring().match;

		if (set_affine_score) {
			mapping.score_type = ScoringType::Affine;
		}
	
		return;
	}

	LevMyers* levMyers;
	levMyers = (sequence_len < 256)
		? (sequence_len < 128) ? levMyers128 : levMyers256
		: levMyers64;

/*	uint32_t flank = (params->sensitive_mode ? params->max_no_errors * params->sensitivity_factor : params->max_no_errors) + 1;
	uint32_t start_pos = ((mapping.dir == genome_t::direct) ? mapping.raw_pos : mapping.raw_pos - sequence_len) - flank;
	uint32_t length = sequence_len + 2 * flank;*/

	uint32_t flank = (params->sensitive_mode ? params->max_no_errors * params->sensitivity_factor : params->max_no_errors) + 1;
	uint32_t start_pos = ((mapping.dir == genome_t::direct) ? mapping.raw_pos : mapping.raw_pos - sequence_len + 1) - flank;
	uint32_t length = sequence_len + 2 * flank + 1;
	uint32_t edit_distance;
	ref_pos_t relative_end_pos;

	levMyers->preprocessRawSeq(sequence, sequence_len, mapping.dir);
	bool ok = levMyers->dynamicProgramming(start_pos, length, mapping.err_edit_distance, relative_end_pos, edit_distance);


#ifdef LEV_MYERS_ASSERTIONS
	uint32_t pos2, edit_distance2;
	levMyers64->preprocessRawSeq(sequence[r], sequence_len[r], mapping.dir);
	levMyers64->dynamicProgramming(start_pos, length, dist, pos2, edit_distance2);
	ASSERT((pos == pos2) && edit_distance == edit_distance2, "CSamGenerator::store_mapped_pair_reads() DP error (mappings > 0)!");
#endif
	if (ok) {
		if (mapping.dir == genome_t::direct) {
			copy_direct(tmp_read_sequence + 1, sequence, sequence_len);
		}
		else {
			convert_to_rev_comp(tmp_read_sequence + 1, sequence, sequence_len);
		}
		ref_copy_direct(tmp_ref_sequence, reference->GetData(), start_pos, length);

		mp_ext_cigar->Reserve(mapping.ext_cigar);

		double affine_score;
		uint32_t num_events;

		ref_pos_t relative_begin_pos =
			levMyers->getExtCigar(mapping.ext_cigar, tmp_ref_sequence, tmp_read_sequence, quality, 
				relative_end_pos, edit_distance, sequence_len, mapping.dir, insertSizeModel->getScoring(), affine_score, num_events);

		int64_t delta = (int64_t)(start_pos + relative_begin_pos) - ((mapping.dir == genome_t::direct) ? mapping.raw_pos : mapping.raw_pos - sequence_len + 1);
		mapping.raw_pos += delta;
		mapping.ref_seq_pos += delta;
		mapping.ref_length = relative_end_pos - relative_begin_pos;
		mapping.method = distant ? MatchingMethod::FromMappingDistant : MatchingMethod::FromMapping;
		mapping.num_events = num_events;

		if (set_affine_score) {
			mapping.score = affine_score;
			mapping.score_type = ScoringType::Affine;
		}


#ifdef LEV_MYERS_ASSERTIONS
		string q1((char*)mapping.ext_cigar);
		uchar_t ext_cigar2[1001];
		levMyers64->getExtCigar(ext_cigar2, tmp_ref_sequence, tmp_read_sequence, pos, edit_distance, sequence_len[r]);
		string q2((char*)ext_cigar2);
		ASSERT(q1 == q2, "CSamGenerator::store_mapped_pair_reads() CIGAR extraction error (mappings > 0)");
#endif
	}
	else {
		++stat_mapped_pe_errors;
		cerr << "##" << id << "(" << mapping.err_edit_distance << " vs " << edit_distance << ")" << endl;
	}
}

// EOF
