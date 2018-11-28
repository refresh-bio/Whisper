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


#ifndef _TYPES_H
#define _TYPES_H

#include "../common/utils.h"
#include <vector>
#include <utility>
#include <limits>
#include <functional>

#undef min
#undef max

// ************************************************************************************
//
// ************************************************************************************
struct file_name_no_t
{
	uint32_t file_no;
	string file_name1;
	string file_name2;

	file_name_no_t() :
		file_no(0), file_name1(""), file_name2("")
	{};
	file_name_no_t(uint32_t _file_no, string _file_name1, string _file_name2) :
		file_no(_file_no), file_name1(_file_name1), file_name2(_file_name2)
	{};
};

// ************************************************************************************
//
// ************************************************************************************
struct res_group_t
{
	uint32_t group_id;
	uchar_t* data;
	uint64_t size;

	res_group_t() :
		group_id(), data(nullptr), size(0)
	{};
	res_group_t(uint32_t _group_id, uchar_t* _data, uint64_t _size) :
		group_id(_group_id), data(_data), size (_size)
	{};
};

// ************************************************************************************
//
// ************************************************************************************
struct reads_bin_t
{
	uint32_t bin_num;
	uint32_t bin_id;
	uchar_t* data;
	uint64_t raw_size;
	uint64_t size;
	uint64_t count;
	uint64_t req_memory;

	reads_bin_t() :
		bin_num(), bin_id(), data(nullptr), raw_size(0), size(0), count(0)
	{};
	reads_bin_t(uint32_t _bin_num, uint32_t _bin_id, uchar_t* _data, uint64_t _raw_size, uint64_t _size, uint64_t _count) :
		bin_num(_bin_num), bin_id(_bin_id), data(_data), raw_size(_raw_size), size (_size), count(_count)
	{};
};

// ************************************************************************************
struct fastq_block_t
{
	read_id_t id_range;
	uchar_t* data;
	uint32_t size;

	fastq_block_t() :
		id_range(empty_read_id), data(nullptr), size(0)
	{};

	fastq_block_t(read_id_t _id_range, uchar_t* _data, uint32_t _size) :
		id_range(_id_range), data(_data), size (_size)
	{};
};

// ************************************************************************************
// Data from bin with mapping results
struct results_bin_t
{
	uchar_t *data;
	uint64_t size;

	results_bin_t():
		data(nullptr), size(0)
	{};

	results_bin_t(uchar_t* _data, uint64_t _size):
		data(_data), size(_size)
	{};
};

// ************************************************************************************
// Whole mapping results for single bin and a few FASTQ blocks
struct mapped_reads_t
{
	results_bin_t results;
	vector<fastq_block_t> fastq_blocks;
	bool all_blocks;
	bool single_end;

	mapped_reads_t() :
		all_blocks(false), single_end(true)
	{};
};

// ************************************************************************************
// Parts of SAM file
struct sam_block_t
{
	uchar_t *data;
	uint64_t size;
	sam_results_t type;

	sam_block_t():
		data(nullptr), size(0), type(sam_results_t::undefined)
	{};

	sam_block_t(uchar_t* _data, uint64_t _size, sam_results_t _type):
		data(_data), size(_size), type(_type)
	{};
};

// ************************************************************************************
// Description of sizes of joined parts of reference sequence
struct seq_desc_t {
	string name;						// sequence name
	int32_t id;							// identifier
	uint64_t size;						// total length of the sequence
	uint64_t pos_in_ref_seq;			// position of the sequence in the compiled reference sequence
	uint64_t no_initial_Ns;				// no. of removed initial Ns in the original sequence
	uint64_t no_final_Ns;				// no. of removed final Ns in the original sequence
	uint64_t no_starting_Ns;			// no. of inserted starting Ns in the compiled reference sequence

	seq_desc_t(const string& _name, int32_t _id, const uint64_t _size, const uint64_t _pos_in_ref_seq,
		const uint64_t _no_initial_Ns, const uint64_t _no_final_Ns, const uint64_t _no_starting_Ns) :
	name(_name), id(_id), size(_size), pos_in_ref_seq(_pos_in_ref_seq), no_initial_Ns(_no_initial_Ns), no_final_Ns(_no_final_Ns), no_starting_Ns(_no_starting_Ns) 
	{};
};

// ************************************************************************************
// Match decription type
class MatchingMethod {
public:
	enum Enumeration { Perfect, PerfectDistant, FromMapping, FromMappingDistant, LevMyers, Clipping, Unmatched };
	
	MatchingMethod() {}
	MatchingMethod(Enumeration v) : value(v) {}
	
	static size_t size() { return 7; }

	std::string toString() {
		switch (value) {
			case Perfect: return "Perfect";
			case PerfectDistant: return "PerfectDistant";
			case FromMapping: return "FromMapping";
			case FromMappingDistant: return "FromMappingDistant";
			case LevMyers: return "LevMyers";
			case Clipping: return "Clipping";
			case Unmatched: return "Unmatched";
			}
		return "invalid";
	}

	bool operator==(MatchingMethod::Enumeration ref) { return value == ref; }
	bool operator!=(MatchingMethod::Enumeration ref) { return value != ref; }
	MatchingMethod& operator=(MatchingMethod::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};

enum class ScoringType : uint16_t { Linear, Affine };

// ************************************************************************************
// Mapping description type
struct mapping_desc_t 
{
	ref_pos_t raw_pos;
	int32_t ref_seq_id;
	int32_t ref_seq_pos;
	genome_t dir;
	uint16_t read_length;
	uint16_t ref_length;
	uint16_t err_edit_distance;
	uint16_t num_events;
	uint16_t left_clipping;
	uint16_t right_clipping;
	ScoringType score_type;
	uchar_t *ext_cigar;
	MatchingMethod::Enumeration method;
	double score;
	uint8_t mapq;
	
	mapping_desc_t() {}

	mapping_desc_t(
		ref_pos_t _raw_pos,
		int32_t _ref_seq_id,
		int32_t _ref_seq_pos,
		genome_t _dir,
		uint16_t _read_length,
		uint16_t _err_edit_distance,
		uchar_t* _ext_cigar) :

		raw_pos(_raw_pos), ref_seq_id(_ref_seq_id), ref_seq_pos(_ref_seq_pos), dir(_dir),
		read_length(_read_length), ref_length(-1), err_edit_distance(_err_edit_distance), num_events(_err_edit_distance),
		left_clipping(0), right_clipping(0), 
		score_type(ScoringType::Linear),
		ext_cigar(_ext_cigar),
		method(MatchingMethod::Unmatched), 
		//linear_score(-std::numeric_limits<double>::max()),
		//affine_score(-std::numeric_limits<double>::max()),
		score(-1),
		mapq(255)
	{};

	static bool compareByScoresDescending(const mapping_desc_t& p, const mapping_desc_t& q) {
		ASSERT(p.score_type == q.score_type, "mapping_desc_t::compareByScoresDescending");
		return (p.score == q.score) ? compareByPositions(p, q) : p.score > q.score;
	}

	static bool compareByScoresDescending_ptr(const mapping_desc_t* p, const mapping_desc_t* q) {
		ASSERT(p->score_type == q->score_type, "mapping_desc_t::compareByScoresDescending_ptr");
		
		return (p->score == q->score) ? compareByPositions_ptr(p, q) : p->score > q->score;
	}
	
	static bool compareByPositions(const mapping_desc_t& p, const mapping_desc_t& q) {
		return (p.ref_seq_id == q.ref_seq_id)
			? p.ref_seq_pos < q.ref_seq_pos
			: p.ref_seq_id < q.ref_seq_id;
	};

	static bool compareByPositions_ptr(const mapping_desc_t* p, const mapping_desc_t* q) {
		return (p->ref_seq_id == q->ref_seq_id)
			? p->ref_seq_pos < q->ref_seq_pos
			: p->ref_seq_id < q->ref_seq_id;
	};

	static bool equalPositions(const mapping_desc_t& p, const mapping_desc_t& q) {
		return p.ref_seq_id == q.ref_seq_id && p.ref_seq_pos == q.ref_seq_pos;
	}
};

// ************************************************************************************
// Mapping pair description type
class mapping_pair_t : public std::pair<mapping_desc_t*, mapping_desc_t*> 
{
public:
	static const int MAX_INSERT_SIZE = std::numeric_limits<int>::max();;

	double score;
	
	mapping_pair_t() : std::pair<mapping_desc_t*, mapping_desc_t*>(), score(-1) {}
	mapping_pair_t(mapping_desc_t* first, mapping_desc_t* second) : std::pair<mapping_desc_t*, mapping_desc_t*>(first, second), score(-1) {}

	mapping_desc_t& operator [](int i) { return i == 0 ? *first : *second; }
	const mapping_desc_t& operator [](int i) const { return i == 0 ? *first : *second; }

	static int calculateInsertSize(const mapping_desc_t& first, const mapping_desc_t& second) {
		// set to infinity if hits are from different chromosomes
		if (first.ref_seq_id != second.ref_seq_id) {
			return MAX_INSERT_SIZE;
		}

		const mapping_desc_t* left = nullptr;
		const mapping_desc_t* right = nullptr;
		
		if (first.ref_seq_pos < second.ref_seq_pos) {
			left = &first;
			right = &second;
		} else if (first.ref_seq_pos == second.ref_seq_pos) {
			if (first.dir == genome_t::direct) { // if positions are equal - treat direc as left
				left = &first;
				right = &second;
			}
			else {
				left = &second;
				right = &first;
			}
		}
		else {
			left = &second;
			right = &first;
		}

		// support reads with proper orientation only
		if (first.dir == second.dir || left->dir == genome_t::rev_comp) {
			return MAX_INSERT_SIZE;
		}

		int out = right->ref_seq_pos - left->ref_seq_pos + right->ref_length;
		return out;
	}

	int calculateInsertSize() const {
		return calculateInsertSize(*first, *second);
	}

	static bool compareByPositions(const mapping_pair_t& p, const mapping_pair_t& q) {
		
		if (p.first->ref_seq_id == q.first->ref_seq_id) {				// same first chromosome
			if (p.first->ref_seq_pos == q.first->ref_seq_pos) {			// same first position	
				if (p.second->ref_seq_id == q.second->ref_seq_id) {		// same second chromosome
					return p.second->ref_seq_pos < q.second->ref_seq_pos;
				}
				return p.second->ref_seq_id < q.second->ref_seq_id;
			}
			return p.first->ref_seq_pos < q.first->ref_seq_pos;
		}
		
		return  p.first->ref_seq_id < q.first->ref_seq_id; 
	};


	static bool compareByScoresDescending(const mapping_pair_t& p, const mapping_pair_t& q) {
		ASSERT(
			p.first->score_type == p.second->score_type && 
			q.first->score_type == q.second->score_type && 
			p.first->score_type == q.first->score_type,
			"mapping_pair_t::compareByScoresDescending");

		return (p.score == q.score)
			? compareByPositions(p, q)
			: p.score > q.score;
	}

	static bool equalPositions(const mapping_pair_t& p, const mapping_pair_t& q) {
		return 
			p.first->ref_seq_id == q.first->ref_seq_id && p.second->ref_seq_id == q.second->ref_seq_id &&
			p.first->ref_seq_pos == q.first->ref_seq_pos && p.second->ref_seq_pos == q.second->ref_seq_pos;
	}
};

// ************************************************************************************
// Evolution model
class scoring_t 
{
public:
	double match;
	
	// linear model - common score for mismatch, gap open, and gap extension
	double error;

	// affine model
	double mismatch;
	double gap_open;
	double gap_extend;

	double clipping;

	scoring_t(double match, double error, double mismatch, double gap_open, double gap_extend, double clipping) :
		match(match), error(error), mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend), clipping(clipping) {}


	const double getOptimisticScore(int read_length, int edit_distance) const {
		double score = (read_length) * match;

		score += std::max(edit_distance * mismatch, gap_open + (edit_distance - 1) * gap_extend);

		return score;
	}
};

// ************************************************************************************
// String read type. Allows reading null terminated strings forward or backward
class string_reader_t
{
public:
	string_reader_t(const uchar_t* s, int len, genome_t dir) : s(s), len(len) {
		if (dir == genome_t::direct) {
			accessor = [this](int i)->uchar_t { return this->s[i];  };
		}
		else {
			accessor = [this](int i)->uchar_t { return this->s[this->len - 1 - i];  };
		}
	}

	uchar_t operator[](int i) { return accessor(i);  }

private:
	const uchar_t *s;
	int len;
	std::function<uchar_t(int)> accessor;
};

#endif

// EOF
