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
		bin_num(0), bin_id(0), data(nullptr), raw_size(0), size(0), count(0), req_memory(0)
	{};
	reads_bin_t(uint32_t _bin_num, uint32_t _bin_id, uchar_t* _data, uint64_t _raw_size, uint64_t _size, uint64_t _count) :
		bin_num(_bin_num), bin_id(_bin_id), data(_data), raw_size(_raw_size), size (_size), count(_count), req_memory(0)
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

	seq_desc_t() :
		name(""), id(-1), size(0), pos_in_ref_seq(0), no_initial_Ns(0), no_final_Ns(0), no_starting_Ns(0)
	{};

	seq_desc_t(const string& _name, int32_t _id, const uint64_t _size, const uint64_t _pos_in_ref_seq,
		const uint64_t _no_initial_Ns, const uint64_t _no_final_Ns, const uint64_t _no_starting_Ns) :
	name(_name), id(_id), size(_size), pos_in_ref_seq(_pos_in_ref_seq), no_initial_Ns(_no_initial_Ns), no_final_Ns(_no_final_Ns), no_starting_Ns(_no_starting_Ns) 
	{};
};

// ************************************************************************************
// Match decription type
class MatchingMethod {
public:
	enum class Enumeration { Perfect, PerfectDistant, FromMapping, FromMappingDistant, LevMyers, Clipping, ClippingInsert, ClippingDelete, 
		VariantInsertLong, VariantDeleteLong, VariantInsertShort, VariantDeleteShort, VariantIndelSV, 
		MappingInsert, MappingDelete, MappingInsertDelete, MappingInsertInsert, MappingDeleteInsert, MappingDeleteDelete,
		MappingClippingMismatches, MappingMismatchesClipping,
		MappingClippingInsert, MappingInsertClipping,
		MappingClippingDelete, MappingDeleteClipping,
		MappingClippingClipping,
		Unmatched };
	
	MatchingMethod() : value(Enumeration::Unmatched) {}
	MatchingMethod(Enumeration v) : value(v) {}
	
	static size_t size() { return 27; }

	std::string toString() {
		switch (value) {
			case Enumeration::Perfect: return "Perfect";
			case Enumeration::PerfectDistant: return "PerfectDistant";
			case Enumeration::FromMapping: return "FromMapping";
			case Enumeration::FromMappingDistant: return "FromMappingDistant";
			case Enumeration::LevMyers: return "LevMyers";
			case Enumeration::Clipping: return "Clipping";
			case Enumeration::ClippingInsert: return "ClippingInsert";
			case Enumeration::ClippingDelete: return "ClippingDelete";
			case Enumeration::VariantInsertLong: return "VariantInsertLong";
			case Enumeration::VariantDeleteLong: return "VariantDeleteLong";
			case Enumeration::VariantInsertShort: return "VariantInsertShort";
			case Enumeration::VariantDeleteShort: return "VariantDeleteShort";
			case Enumeration::VariantIndelSV: return "VariantIndelSV";
			case Enumeration::MappingInsert: return "MappingInsert";
			case Enumeration::MappingDelete: return "MappingDelete";
			case Enumeration::MappingInsertDelete: return "MappingInsertDelete";
			case Enumeration::MappingInsertInsert: return "MappingInsertInsert";
			case Enumeration::MappingDeleteInsert: return "MappingDeleteInsert";
			case Enumeration::MappingDeleteDelete: return "MappingDeleteDelete";
			case Enumeration::MappingMismatchesClipping: return "MappingMismatchesClipping";
			case Enumeration::MappingClippingMismatches: return "MappingClippingMismatches";
			case Enumeration::MappingInsertClipping: return "MappingInsertClipping";
			case Enumeration::MappingClippingInsert: return "MappingClippingInsert";
			case Enumeration::MappingDeleteClipping: return "MappingDeleteClipping";
			case Enumeration::MappingClippingDelete: return "MappingClippingDelete";
			case Enumeration::MappingClippingClipping: return "MappingClippingClipping";
			case Enumeration::Unmatched: return "Unmatched";
			}
		return "invalid";
	}

	bool operator==(MatchingMethod::Enumeration ref) { return value == ref; }
	bool operator!=(MatchingMethod::Enumeration ref) { return value != ref; }
	MatchingMethod& operator=(MatchingMethod::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};

// ************************************************************************************
struct candidate_mapping_t
{
	static const uint32_t max_size = 11;
	static const uint32_t min_size = 5;

	mapping_type_t type;
	ref_pos_t pos;
	genome_t direction;
	uint32_t no_mismatches;
	uint32_t matching_len1;
	uint32_t matching_len2;
	int32_t indel_len1;
	int32_t indel_len2;
	uint32_t clipping_len;

	double penalty;

private:
	uint32_t calc_penalty_lev() { return no_mismatches; }
	uint32_t calc_penalty_mismatches_clipping() { return no_mismatches + 1; }
	uint32_t calc_penalty_indel_clipping() { return no_mismatches + 2; }
	uint32_t calc_penalty_indel1() { return no_mismatches + 1; }
	uint32_t calc_penalty_indel2() { return no_mismatches + 2; }
	uint32_t calc_penalty_clipping_clipping() { return no_mismatches + 2; }

public:
	candidate_mapping_t() :
		type(mapping_type_t::none),
		pos(0),
		direction(genome_t::direct),
		no_mismatches(0),
		matching_len1(0),
		matching_len2(0),
		indel_len1(0),
		indel_len2(0),
		clipping_len(0),
		penalty(0)
	{};

	candidate_mapping_t(const candidate_mapping_t& mapping) = default;

	bool operator<(const candidate_mapping_t& mapping);
	candidate_mapping_t& operator=(const candidate_mapping_t& mapping);

	// type: lev
	template<mapping_type_t mapping_type>
	static candidate_mapping_t construct(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches);

	// type: mismatches_clipping, clipping_mismatches
	template<mapping_type_t mapping_type>
	static candidate_mapping_t construct(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, 
		uint32_t _clipping_len);

	// type: indel_clipping, clipping_indel
	template<mapping_type_t mapping_type>
	static candidate_mapping_t construct(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, 
		uint32_t _matching_len1, int32_t _indel_len1, uint32_t _matching_len2);

	// type: indel1
	template<mapping_type_t mapping_type>
	static candidate_mapping_t construct(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, 
		uint32_t _matching_len1, int32_t _indel_len1);

	// type: indel2
	template<mapping_type_t mapping_type>
	static candidate_mapping_t construct(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, 
		uint32_t _matching_len1, int32_t _indel_len1, uint32_t _matching_len2, int32_t _indel_len2);

	void construct(candidate_mapping_t& mapping);

	void calc_penalty();

	uchar_t indel_size_encode(int32_t x) const
	{
		if (x > 0)
			return (uchar_t)x;
		else
			return (uchar_t)(-x + 128);
	}

	int32_t indel_size_decode(uchar_t x) const
	{
		if (x < 128)
			return (int32_t)x;
		else
			return -(int32_t)(x - 128);
	}

	size_t serialize(uchar_t* ptr);
	size_t deserialize(uchar_t* ptr);
	static size_t check_rec_size(uchar_t* ptr);
	static uchar_t* skip(uchar_t* ptr);

	static candidate_mapping_t better(candidate_mapping_t& x, candidate_mapping_t& y);
};


enum class ScoringType : uint16_t { Linear, Affine };

// ************************************************************************************
// Mapping description type
struct mapping_desc_t 
{
	candidate_mapping_t mapping;
	int32_t ref_seq_id;
	int32_t ref_seq_pos;
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
	
	uint16_t indel_to_mismatch_score(int32_t size, double gap_ins_open, double gap_ins_extend, double gap_del_open, double gap_del_extend, double mismatch_score)
	{
		double affine_score = 0;

		if (size > 0)
			affine_score = gap_del_open + ((double)size - 1) * gap_del_extend;
		else
			affine_score = gap_ins_open + ((double)-size - 1) * gap_ins_extend;

		return (uint16_t)(affine_score / mismatch_score);
	}

	mapping_desc_t() :
		mapping(),
		ref_seq_id(0),
		ref_seq_pos(0),
		read_length(0),
		ref_length(0),
		err_edit_distance(0),
		num_events(0),
		left_clipping(0),
		right_clipping(0),
		score_type(ScoringType::Linear),
		ext_cigar(nullptr),
		method(MatchingMethod::Enumeration::Unmatched),
		score(0.0),
		mapq(0)
	{}

	mapping_desc_t(
		candidate_mapping_t _mapping,
		int32_t _ref_seq_id,
		int32_t _ref_seq_pos,
		uint16_t _read_length,
		uchar_t* _ext_cigar,
		double gap_ins_open, double gap_ins_extend, double gap_del_open, double gap_del_extend, double mismatch_score) :
		mapping(_mapping), ref_seq_id(_ref_seq_id), ref_seq_pos(_ref_seq_pos), 
		read_length(_read_length), ref_length(-1), 
		err_edit_distance(0), num_events(0),
		left_clipping(0), right_clipping(0),
		score_type(ScoringType::Linear),
		ext_cigar(_ext_cigar),
		method(MatchingMethod::Enumeration::Unmatched), 
		//linear_score(-std::numeric_limits<double>::max()),
		//affine_score(-std::numeric_limits<double>::max()),
		score(-1),
		mapq(255)
	{
		switch (mapping.type)
		{
		case mapping_type_t::lev:
			err_edit_distance = mapping.no_mismatches;
			num_events = mapping.no_mismatches;
			break;
		case mapping_type_t::indel1:
			err_edit_distance = mapping.no_mismatches +
				indel_to_mismatch_score(mapping.indel_len1, gap_ins_open, gap_ins_extend, gap_del_open, gap_del_extend, mismatch_score);
			num_events = mapping.no_mismatches + 1;
			break;
		case mapping_type_t::indel2:
			err_edit_distance = mapping.no_mismatches +
				indel_to_mismatch_score(mapping.indel_len1, gap_ins_open, gap_ins_extend, gap_del_open, gap_del_extend, mismatch_score) +
				indel_to_mismatch_score(mapping.indel_len2, gap_ins_open, gap_ins_extend, gap_del_open, gap_del_extend, mismatch_score);
			num_events = mapping.no_mismatches + 2;
			break;
		case mapping_type_t::clipping_clipping:
			err_edit_distance = mapping.no_mismatches + 2;
			num_events = mapping.no_mismatches + 2;
			break;
		case mapping_type_t::clipping_indel:
		case mapping_type_t::indel_clipping:
			err_edit_distance = mapping.no_mismatches + 1 +
				indel_to_mismatch_score(mapping.indel_len1, gap_ins_open, gap_ins_extend, gap_del_open, gap_del_extend, mismatch_score);
			num_events = mapping.no_mismatches + 2;
			break;
		case mapping_type_t::mismatches_clipping:
		case mapping_type_t::clipping_mismatches:
			err_edit_distance = mapping.no_mismatches + 1;
			num_events = mapping.no_mismatches + 1;
			break;
		}
	};

	static bool compareByScoresDescending(const mapping_desc_t& p, const mapping_desc_t& q) {
		ASSERT(p.score_type == q.score_type, "mapping_desc_t::compareByScoresDescending");
		return p.score > q.score;
	}

	static bool compareByScoresAndPosDescending(const mapping_desc_t& p, const mapping_desc_t& q) {
		ASSERT(p.score_type == q.score_type, "mapping_desc_t::compareByScoresDescending");
		return (p.score == q.score) ? compareByPositions(p, q) : p.score > q.score;
	}

	static bool compareByScoresDescending_ptr(const mapping_desc_t* p, const mapping_desc_t* q) {
		ASSERT(p->score_type == q->score_type, "mapping_desc_t::compareByScoresDescending_ptr");
		
		return p->score > q->score;
	}
	
	static bool compareByScoresAndPosDescending_ptr(const mapping_desc_t* p, const mapping_desc_t* q) {
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
	bool is_distant;
	
	mapping_pair_t() : std::pair<mapping_desc_t*, mapping_desc_t*>(), score(-1), is_distant(false) {}
	mapping_pair_t(mapping_desc_t* first, mapping_desc_t* second, bool _is_distant = false) : std::pair<mapping_desc_t*, mapping_desc_t*>(first, second), score(-1), is_distant(_is_distant) {}

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
			if (first.mapping.direction == genome_t::direct) { // if positions are equal - treat direct as left
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
		if (first.mapping.direction == second.mapping.direction || left->mapping.direction == genome_t::rev_comp) {
			return MAX_INSERT_SIZE;
		}

		return right->ref_seq_pos - left->ref_seq_pos + right->ref_length;
	}

	int calculateInsertSize() const {
		return calculateInsertSize(*first, *second);
	}

	static bool compareByPositions(const mapping_pair_t& p, const mapping_pair_t& q) {		
		if (p.first->ref_seq_id == q.first->ref_seq_id) {				// same first chromosome
			if (p.first->ref_seq_pos == q.first->ref_seq_pos) {			// same first position	
				if (p.second->ref_seq_id == q.second->ref_seq_id) {		// same second chromosome
//					return p.second->ref_seq_pos < q.second->ref_seq_pos;
					if(p.second->ref_seq_pos != q.second->ref_seq_pos)
						return p.second->ref_seq_pos < q.second->ref_seq_pos;
					else {
						if (p.is_distant == q.is_distant)
							return false;
						else
							return q.is_distant;
					}						
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

		return p.score > q.score;
	}

	static bool compareByScoresAndPosDescending(const mapping_pair_t& p, const mapping_pair_t& q) {
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
//	double gap_open;
//	double gap_extend;
	double gap_ins_open;
	double gap_ins_extend;
	double gap_del_open;
	double gap_del_extend;

	double clipping;

	scoring_t(double match, double error, double mismatch, 
//		double gap_open, double gap_extend, 
		double gap_ins_open, double gap_ins_extend, double gap_del_open, double gap_del_extend,
		double clipping) :
		match(match), error(error), mismatch(mismatch), 
//		gap_open(gap_open), gap_extend(gap_extend), 
		gap_ins_open(gap_ins_open), gap_ins_extend(gap_ins_extend), gap_del_open(gap_del_open), gap_del_extend(gap_del_extend),
		clipping(clipping) {}


	const double getOptimisticScore(int read_length, int edit_distance) const {
		double score = (read_length) * match;

		score += std::max(edit_distance * mismatch, 
			std::max(gap_ins_open + (edit_distance - 1.0) * gap_ins_extend, 
				gap_del_open + (edit_distance - 1.0) * gap_del_extend));

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

// ************************************************************************************
//typedef enum class variant_type {snp, short_indel, long_indel};

// ************************************************************************************
// Description of referential variants
struct variant_snp_t
{
	uint8_t is_ref;
	uint8_t symbol;
	uint16_t freq;
	uint32_t pos;

	variant_snp_t() :
		is_ref(0), symbol(0), freq(0), pos(0)
	{};

	variant_snp_t(uint8_t _is_ref, uint8_t _symbol, uint16_t _freq, uint32_t _pos) :
		is_ref(_is_ref), symbol(_symbol), freq(_freq), pos(_pos)
	{};
};

// ************************************************************************************
struct variant_ins_t
{
	uint32_t ins_desc_pos;
	uint32_t ins_len;
	uint16_t freq;
	uint32_t pos;

	variant_ins_t() :
		ins_desc_pos(0), ins_len(0), freq(0), pos(0)
	{};

	variant_ins_t(uint32_t _ins_desc_pos, uint32_t _ins_len, uint16_t _freq, uint32_t _pos) :
		ins_desc_pos(_ins_desc_pos), ins_len(_ins_len), freq(_freq), pos(_pos)
	{};
};

// ************************************************************************************
struct variant_del_t
{
	uint32_t del_len;
	uint16_t freq;
	uint32_t pos;

	variant_del_t() :
		del_len(0), freq(0), pos(0)
	{};

	variant_del_t(uint32_t _del_len, uint16_t _freq, uint32_t _pos) :
		del_len(_del_len), freq(_freq), pos(_pos)
	{};
};

// ************************************************************************************
struct variant_sv_t
{
	uint32_t ins_desc_pos;
	uint32_t ins_len;
	uint32_t del_len;
	uint16_t freq;
	uint32_t pos;

	variant_sv_t() :
		ins_desc_pos(0), ins_len(0), del_len(0), freq(0), pos(0)
	{};

	variant_sv_t(uint32_t _ins_desc_pos, uint32_t _ins_len, uint32_t _del_len, uint16_t _freq, uint32_t _pos) :
		ins_desc_pos(_ins_desc_pos), ins_len(_ins_len), del_len(_del_len), freq(_freq), pos(_pos)
	{};
};

#endif

// EOF
