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


#ifndef _MAPPING_CORE_H
#define _MAPPING_CORE_H

#include "../common/defs.h"
#include "../common/idstore.h"
#include "../common/queue.h"
#include "../common/mmgr.h"
#include "../common/stats.h"
#include "../common/timer.h"
#include "../common/results.h"
#include "../common/params.h"
#include "edit_dist.h"
#include "reads.h"
#include "../libs/vectorclass.h"
#include "LevMyers.h"
#include <vector>

#define SELECT_NON_ZERO_VALUES
#define USE_RESULT_OF_A_MEMCMP

#define USE_128b_SSE_CF_AND_DP

#define SWAP_NIBBLES(x)		((((x) & 0xf) << 4) + ((x) >> 4))
#define HI_NIBBLE(x)		((x) >> 4)
#define LO_NIBBLE(x)		((x) & 0x0f)
#define RAW_HI_NIBBLE(x)	((x) & 0xf0)

#define CMP_DIFFERENT_SWAPPED_NIBBLES(x, y)	((LO_NIBBLE(x) < LO_NIBBLE(y)) ? 1 :\
											 (LO_NIBBLE(x) > LO_NIBBLE(y)) ? -1 :\
											 ((x) < (y)) ? 1 : -1)

using namespace std;

struct CAux_arrays {
	std::vector<uint32_t> howManySubstr_array;
	std::vector<uint32_t> whereToPut_array;
	std::vector<uint32_t> localization_array;
};

struct CAux_params {
	uint32_t substring_length;
	vector<uint32_t> *substr_pos_arr;
	vector<int32_t> *substr_pos_in_text_arr;
	vector<pair<uint32_t, uint32_t>> *substr_pos_mapping_arr;
	uint32_t LEN;
	uint32_t no_of_substrings_for_DP;
	uint32_t no_of_substrings_in_text_for_DP;
};

class CMappingCore
{
	int64_t no_of_CF_accepted;
	int64_t no_of_CF_discarded;
	int64_t no_of_CF_accepted_in_malicious_group;
	int64_t no_of_CF_discarded_in_malicious_group;
	int64_t no_of_Lev_positive;

	static uchar_t cmp_lut[256];
	static uchar_t rc_lut[80];
	static uint32_t cf_lut[256];
	static uint64_t CF_lut_64b[256][2];
	static uchar_t list[17][16];
	static uint32_t CF_lut_for_SSE[256];
	static uint32_t dense_symbols[256];
	Vec16c CF_SSE_lut[256];

	vector<uint32_t> prefix_map;
	CRegisteringQueue<reads_bin_t> *q_bins_read;
	CRegisteringQueue<reads_bin_t> *q_bins_write;
	CRegisteringQueue<res_group_t> *q_map_res;
	CMemoryMonitor *mem_monitor;
	CMemoryPool<uchar_t> *mp_bins_write;
	CMemoryPool<uchar_t> *mp_map_res;
	CRunningStats *running_stats;
	CStopWatch watch;
	uint32_t stage_id;
	uint32_t stage_major, stage_minor;
	mapping_mode_t mapping_mode;

	CReadsCollector *reads_collector;
	CReadsDeliverer *reads_deliverer;
	CMappingResultsCollector *results_collector;

	CReference *reference;
	CSuffixArray *sa_dir;
	CSuffixArray *sa_rc;
	bool sa_in_ram;

	CProgress *progress;

	CEditDist *edit_dist;
	LevMyers* levMyers64;
	LevMyers* levMyers128;
	LevMyers* levMyers256;

	uint32_t no_bins;
	uint32_t bin_prefix;
	uint32_t read_id_len;
	uint32_t start_pos;
	uint32_t end_pos;
	uint32_t no_groups;
	uint32_t verbosity_level;
	uint32_t next_start_pos;
	uint32_t next_end_pos;

	uint64_t n_reads;
	uint64_t n_parts;
	uint64_t max_read_len;
	uint32_t min_read_len;

	bool sensitive_mode;
	double sensitivity_factor;
	uint32_t max_CF_difference;
	uint32_t max_no_mappings;

	vector<string> prefixes;
	uint32_t prefix_len;
	uint32_t sa_prefix_overhead;
	uint32_t sa_dir_bin_size;		// size of SA (dir) related to current bin
	uint32_t sa_rc_bin_size;		// size of SA (rc) related to current bin
	uint32_t n_reads_in_bin;		// no. of reads in current bin

	// ---------------------------------
	vector<pair<uint32_t, uint32_t>> segment_segments;
	CParams *params_aux_ptr;
	uchar_t *ref_ptr;
	uint32_t ref_size;
	uint32_t predef_match_length;
	uchar_t* cmp_lut_ptr;
	uint32_t CUR_OFFSET;

	bool dir_inside_malicious_group;
	bool rc_inside_malicious_group;
	bool dir_aux_struct_valid;
	bool rc_aux_struct_valid;
	CAux_arrays* dir_aux_arr;
	CAux_arrays* rc_aux_arr;
	CAux_params* dir_aux_params;
	CAux_params* rc_aux_params;
	uint32_t malicious_group_length;
	uint32_t no_of_substrings;
	uchar_t* old_sequence;
	uchar_t* old_sequence_sft;
	vector<vector<uint32_t>> substr_pos_all_lengths;
	vector<vector<int32_t>> substr_pos_in_text_all_lengths;
	vector<vector<pair<uint32_t, uint32_t>>> substr_pos_mapping_all_lengths;

	vector<uint32_t> group_size_substr_length_mapping;
	vector<pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator>> in_list_desc;

	double sa_dir_frac;
	double sa_rc_frac;
	double sa_frac;
	uint32_t sa_dir_skip;
	uint32_t sa_rc_skip;
	uint32_t* CF_array_dir_gen;
	vector<pair<uint64_t, uint64_t>> CF_vector_64b_dir_gen;
	vector<pair<uint64_t, uint64_t>> CF_vector_64b_rc_gen;

	vector<Vec16c> CF_vector_128b_SSE_dir_gen;
	vector<Vec16c> CF_vector_128b_SSE_rc_gen;

	uint32_t* CF_array_rc_gen;
	uint32_t CF_value_for_a_dir_pattern;
	uint32_t CF_value_for_a_rc_pattern;
	uint64_t CF_value_64b_for_a_dir_pattern[2];
	Vec16c CF_value_SSE_for_a_dir_pattern;
	Vec16c CF_value_SSE_for_a_rc_pattern;
	uint64_t CF_value_64b_for_a_rc_pattern[2];
	uint32_t* prev_dp_ptr;
	uint32_t* curr_dp_ptr;
	uchar_t* genome_prefetch;

	vector<pair<uint32_t, uint32_t>> dir_match_positions;
	vector<pair<uint32_t, uint32_t>> rc_match_positions;

	//	vector<tuple<read_id_t, uint32_t, genome_t, uint32_t>> collected_mappings;
	CMappingsHeapGatherer collected_mappings;

	uint32_t* cur_substr_values_from_pattern;
	vector<uint32_t> sum_of_localizations_dir, sum_of_localizations_rc;
	vector<uchar_t> list_merger_tmp_dir, list_merger_tmp_rc;
	vector<uchar_t> identical_sequences_dir, identical_sequences_rc;
	uint32_t old_pos_in_text_dir, old_pos_in_text_rc;


	//****************************************************************************
	//*********************** mapping_core.cpp ***********************************
	//****************************************************************************
	inline uint32_t get_next_bin_id(uchar_t* data);
	inline uint32_t get_read_prefix(uchar_t* data, uint32_t size);

	void copy_reads(reads_bin_t bin);
	void complete();
	void process_reads_init_stage(reads_bin_t bin);
	void process_reads(reads_bin_t bin);
	void reverse_read(uchar_t *pattern_ptr, uint32_t pattern_len,
		uchar_t **rc_pattern_ptr, uchar_t **rc_pattern_ptr_sft, uchar_t *reverse_lut);

	//****************************************************************************
	//*********************** mapping_core_str_comp.cpp **************************
	//****************************************************************************
	void    string_copy(uchar_t *target, uchar_t *source, uint32_t len);
	int32_t compare_dir_str_with_gen(uchar_t* text_ptr, uchar_t *pattern_ptr, uint32_t len_to_compare, bool odd_pos);
	int32_t compare_rc_str_backward(uchar_t* text_ptr, uchar_t *pattern_ptr, uint32_t len_to_compare, bool odd_pos);

	bool    compare_str_to_str(uchar_t *pattern_ptr, uint32_t pattern_len, uchar_t *old_pattern_ptr);
	bool    compare_str_to_str_in_text(uint32_t new_seq_start_pos, uint32_t len, uint32_t old_seq_start_pos);
	int32_t check_basic_part_of_diff_len_reads(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
		uint32_t position_in_text, bool dir_gen_flag);
	int32_t check_matching_part(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
		uint32_t position_in_text, bool dir_gen_flag);

	//****************************************************************************
	//*********************** mapping_core_substrings.cpp **************************
	//****************************************************************************
	void	 adjust_substr_pos_arr();
	void	 adjust_group_sizes(reads_bin_t bin);
	uint32_t get_substring_value_from_pattern(uint32_t substring_length, uchar_t *pattern, uint32_t k);
	uint32_t get_substring_value(uint32_t substring_length, int32_t start_substring_position, uint32_t ref_index_given_in_chars);
	vector<uint32_t>::iterator list_merger(vector<pair<vector<uint32_t>::iterator, vector<uint32_t>::iterator>> in_lists, uint32_t min_occ,
		vector<uchar_t> *list_merger_tmp, vector<uint32_t>::iterator out_iter);


	//****************************************************************************
	//*********************** mapping_core_dp.cpp ********************************
	//****************************************************************************
	uint32_t search_in_malicious_group_with_DP(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
		uint32_t* sub_SA, uint32_t sub_sa_size,
		CAux_arrays* aux_arr, CAux_params* aux_params, bool* aux_struct_valid, bool dir_gen_flag, uint32_t max_mismatches, uint32_t sub_SA_index_offs);
	void build_aux_structures_for_DP(uint32_t* sub_SA, uint32_t sub_size, uint32_t pattern_len, CAux_arrays* aux_arr, CAux_params* aux_params, bool dir);
	void calculate_CF_128b_SEE_for_dir_and_rc_gen(uint32_t fixedLeft, uint32_t fixedRight, uint32_t pattern_len, uint32_t position_read_from_SA_array, Vec16c& CF, bool dir_gen_flag);
	void calculate_CF_128b_SSE_for_a_pattern(uchar_t *pattern, uchar_t *pattern_sft, uint32_t pattern_len, uint32_t fixedLeft, uint32_t fixedRight, Vec16c& CF_value);
	bool test_CF_128b_SSE_difference(Vec16c CF_1, Vec16c CF_2, uint32_t max_diff);

	//****************************************************************************
	//*********************** mapping_core_sa.cpp ********************************
	//****************************************************************************
	uint32_t search_in_unknown_range_for_exact_match_diff_len_version(uchar_t* pattern_ptr, uchar_t* pattern_ptr_sft, uint32_t pattern_len,
		uint32_t* SA, uint32_t sa_size,
		uint32_t start_SA_index, uint32_t* last_index, bool dir_gen_flag, uint32_t max_mismatches);
	uint32_t search_in_known_range_for_exact_match_diff_len_version(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
		uint32_t* SA, uint32_t sa_size,
		uint32_t start_SA_index, uint32_t last_index, bool dir_gen_flag, uint32_t max_mismatches);
	uint32_t search_in_unknown_range_dir_and_rc(uchar_t* pattern_ptr, uchar_t* pattern_ptr_sft, uint32_t pattern_len,
		uint32_t* SA, uint32_t sa_size,
		uint32_t start_SA_index, uint32_t* last_index, bool dir_gen_flag, uint32_t max_mismatches);
	uint32_t search_in_known_range_dir_and_rc(uchar_t* pattern_ptr, uchar_t* pattern_sft_ptr, uint32_t pattern_len,
		uint32_t* SA, uint32_t sa_size,
		uint32_t start_SA_index, uint32_t last_index, bool dir_gen_flag, uint32_t max_mismatches);
	int32_t binary_sa_search(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len, uint32_t len_to_compare,
		uint32_t* SA, uint32_t sa_size, uint32_t *idx, bool dir_gen_flag);

	//****************************************************************************
	//*********************** mapping_core_diff_len_reads.cpp ********************************
	//****************************************************************************
	bool  findExactMatch_diff_len_version(uchar_t* pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
		bool basic_data_part_is_new,
		uint32_t *sa_last_index, uint32_t* sa_start_index, uint32_t *sa_part, uint32_t sa_size, bool dir_gen_flag);

	uint32_t findMismatches_dir_and_rc(uchar_t *pattern_ptr, uchar_t *pattern_sft_ptr, uint32_t pattern_len,
		bool matching_part_is_new,
		uint32_t *sa_part, uint32_t sa_size,
		uint32_t* sa_dir_last_index, uint32_t* sa_dir_start_index,
		bool dir_gen_flag, bool* inside_malicious_group, bool* aux_struct_valid, uint32_t max_mismatches);
	//--------------------------------->

public:
	CMappingCore(CParams *params, CObjects *objects, uint32_t _stage_id, mapping_mode_t _mapping_mode);
	~CMappingCore();

	void operator()();
};
#endif

// EOF 