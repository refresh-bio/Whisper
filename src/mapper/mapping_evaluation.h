// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.0
// Date    : 2017-12-24
// License : GNU GPL 3
// *******************************************************************************************


#pragma once
#include "../common/types.h"
#include "../common/utils.h"
#include "../mapper/probabilistic.h"

// ************************************************************************************
// Universal evaluation class
template <class T>
struct Evaluator 
{
	double operator()(const T& p) const;
};

// ************************************************************************************
// Smith-Waterman mapping evaluation
template <>
struct Evaluator<mapping_desc_t> 
{
	const InsertSizeModel& model;
	
	Evaluator(const InsertSizeModel& insertSizeModel) : model(insertSizeModel) {}

	double operator()(const mapping_desc_t& p) const {
		return 
			static_cast<float>(p.read_length - p.left_clipping - p.right_clipping - p.err_edit_distance) * model.getScoring().match +
			static_cast<float>(p.err_edit_distance) *  model.getScoring().error;
	}
};

// ************************************************************************************
// Universal pair evaluation
template <>
struct Evaluator<mapping_pair_t> 
{
	// Model of insert size
	const InsertSizeModel& insertSizeModel;

	Evaluator(const InsertSizeModel& insertSizeModel) : insertSizeModel(insertSizeModel) {}

	double operator()(const mapping_pair_t& pair) const {
		int insertSize = pair.calculateInsertSize();
		double s1 = pair.first->score;
		double s2 = pair.second->score;
		double penalty = insertSizeModel.calculatePenalty(insertSize);
		return  s1 + s2 - penalty;		
	}
};

// ************************************************************************************
// Mapping quality
struct MappingEvaluation 
{
	struct result_t {
		double score; 
		uint32_t count;

		static bool compare(const result_t &p, const result_t&q) {
			return p.score > q.score;
		}
	};

	/// Score and count for best result
	result_t best[2];

	double mapq_mult;

	double mapq_div;

	MappingEvaluation(double mapq_mult, double mapq_div) : best{ {0,0}, {0,0} }, mapq_mult(mapq_mult), mapq_div(mapq_div) {}

	MappingEvaluation(const MappingEvaluation& a, const MappingEvaluation& b) {
		if (result_t::compare(a.best[0], b.best[0])) { 
			this->best[0] = a.best[0]; // a[0] > ...

			if (result_t::compare(a.best[1], b.best[0])) {
				this->best[1] = a.best[1]; // a[0] > a[1] > b[0] > ...
			} else if (result_t::compare(b.best[0], a.best[1])) {
				this->best[1] = b.best[0]; // a[0] > b[0] > a[1] > ...
			} else {
				this->best[1] = a.best[1]; // a[0] > a[1] = b[0] > ...
				this->best[1].count += b.best[0].count;
			}
		}
		else if (result_t::compare(b.best[0], a.best[0])) {
			this->best[0] = b.best[0]; // b[0] > ...

			if (result_t::compare(b.best[1], a.best[0])) {
				this->best[1] = b.best[1]; // b[0] > b[1] > a[0] > ...
			}
			else if (result_t::compare(a.best[0], b.best[1])) {
				this->best[1] = a.best[0]; // b[0] > a[0] > b[1] > ...
			}
			else {
				this->best[1] = b.best[1]; // b[0] > b[1] = a[0] > ...
				this->best[1].count += a.best[0].count;
			}
		}
		else {
			this->best[0] = a.best[0]; // a[0] = b[0] > ...
			this->best[0].count = b.best[0].count; 

			if (result_t::compare(a.best[1], b.best[1])) {
				this->best[1] = a.best[1];   // a[0] = b[0] > a[1] > ...
			}
			else if (result_t::compare(b.best[1], a.best[1])) {
				this->best[1] = b.best[1];  // a[0] = b[0] > b[1] > ...
			}
			else {
				this->best[1] = a.best[1]; // a[0] = b[0] > a[1] = b[1]
				this->best[1].count += b.best[1].count;
			}
		}	
	}

	// Calculates MAPQ on the basis of scores
	void initialize(std::vector<mapping_desc_t*>& collection) {
		
		best[0].score = collection[0]->score;
		best[0].count = 1;
		auto prev = collection.begin();
		auto it = std::next(collection.begin());
		while (it != collection.end() && (*it)->score == best[0].score) {
			if (!mapping_desc_t::equalPositions(**it, **prev)) {
				++best[0].count;
			}

			++it;
			++prev;
		}

		if (it != collection.end()) {
			auto it2 = it;
			best[1].score = (*it)->score;
			best[1].count = 1;
			prev = it2;
			++it2;
			while (it2 != collection.end() && (*it2)->score == best[1].score) {
				if (!mapping_desc_t::equalPositions(**it2, **prev)) {
					++best[1].count;
				}

				++it2;
				++prev;
			}
		}

		uint8_t mapq = calculateBestMapQ();

		for (auto jt = collection.begin(); jt != it; ++jt) {
			(*jt)->mapq = mapq;
		}

		for (auto jt = it; jt != collection.end(); ++jt) {
			(*jt)->mapq = 0;
		}
	}

	// Calculates MAPQ on the basis of scores
	void initialize(InsertSizeModel& model, double discretizationThreshold, std::vector<mapping_pair_t>& collection) {
		
		best[0].score = collection[0].score;
		best[0].count = 1;
		auto prev = collection.begin();
		auto it = std::next(collection.begin());
		while (it != collection.end() && it->score > best[0].score - discretizationThreshold) {
			if (!mapping_pair_t::equalPositions(*it, *prev)) {
				++best[0].count;
			}

			++it;
			++prev;
		}

		if (it != collection.end()) {
			best[1].score = it->score;
			best[1].count = 1;
			prev = it;
			++it;
			while (it != collection.end() && it->score > best[1].score - discretizationThreshold) {
				if (!mapping_pair_t::equalPositions(*it, *prev)) {
					++best[1].count;
				}
				
				++it;
				++prev;
			}
		}

		uint8_t pairq = calculateBestMapQ();

		// fixme: hardcoded solution to multiple equally good distant mappings
		auto& top = collection[0];
		if (top.calculateInsertSize() > model.mean() + model.getSaturationDev() && (top.first->mapq <= 3 || top.second->mapq <= 3)) {
			pairq = std::min(top.first->mapq, top.second->mapq);
		}
		
		if (best[0].count == 1) {
			collection[0].first->mapq = std::max(pairq, collection[0].first->mapq);
			collection[0].second->mapq = std::max(pairq, collection[0].second->mapq);
		}
	}

private:
	uint32_t calculateBestMapQ() { 
		if (best[0].count > 1)
			return (best[0].count == 2) ? 3 : (best[0].count == 3) ? 2 : (best[0].count < 10) ? 1 : 0;

		if (!best[1].count)
			return 60;

		int32_t r = 60;

		// error (negative denominator)
		int32_t delta = (int32_t) (mapq_mult * log2(1 + best[1].count) / std::max(best[0].score - best[1].score, mapq_div));
		r -= delta;

		if (r < 0)
			return 0;

		return (uint32_t) r;
	}
};

// EOF
