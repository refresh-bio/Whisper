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


#ifndef _PARAMS_H
#define _PARAMS_H

#include "../common/defs.h"
#include "mmgr.h"
#include "queue.h"
#include "stats.h"
#include "../common/types.h"
#include "../common/reference.h"
#include "sa.h"
#include "../common/variant.h"
#include <string>
#include <vector>
#include <ostream>
#include <sstream>

using namespace std;

// ************************************************************************************
// Command line params
struct CCmdParams
{
	// Files and directories
	string command_line;
	vector<pair<string, string>> input_file_names;
	string index_name;
	string project_name;
	string temp_prefix;
	bool is_fasta;
	bool keep_temporary_files;
	uint32_t gzipped_SAM_level;
	bool store_BAM;
	bool use_stdout;
	string read_group_line;

	bool developer_mode;
	bool log_stats;

	// Memory, disk, threads configuration
	uint32_t no_bins;
	uint32_t no_threads;
	uint64_t max_total_memory;
	bool sa_in_ram;

	// Mapping 
	uint32_t min_read_len;
	uint32_t max_read_len;
	bool paired_reads;
	uint32_t max_no_errors;
	uint32_t max_no_mappings;
	mapping_mode_t mapping_mode;
	mapping_orientation_t mapping_orientation;

	uint32_t mapping_counter_size;
	uint32_t max_mate_distance;
	double max_mate_edit_distance_frac;
	bool sensitive_mode;
	double sensitivity_factor;

	double penalty_saturation_sigmas;
	double clipping_distance_sigmas;
	double high_confidence_sigmas;
	double score_discretization_threshold;
	uint32_t mask_low_quality_bases;
	uint32_t hit_merging_threshold;
	bool hit_merging_wrt_first;
	
	double match_score;
	double mismatch_score;
//	double gap_open;
//	double gap_extend;
	double gap_ins_open;
	double gap_ins_extend;
	double gap_del_open;
	double gap_del_extend;
	double clipping_score;

	uint32_t min_approx_indel_len;
	uint32_t max_approx_indel_len;
	uint32_t max_approx_indel_mismatches;
	double min_clipped_factor;

	double mapq_mult;
	double mapq_div;

	bool enable_boundary_clipping;
	bool enable_paired_clipping;
	bool enable_var_indel_long;
	bool enable_var_indel_short;
	bool enable_var_snp;
	bool enable_short_reads;
	bool enable_mapping_indels;
	bool enable_short_indel_refinement;

	// User interface
	uint32_t verbosity_level;
	string filter;

	CCmdParams()
	{
		developer_mode = false;
		log_stats = false;

		no_bins = 384;
		no_threads = 0;
		max_total_memory = ((uint64_t) 16) << 30;		// 16 GB

		// Mapping 
		max_no_errors = 0;
		mapping_mode = mapping_mode_t::first;
		mapping_orientation = mapping_orientation_t::forward_reverse;
		max_mate_distance = 1000;
		max_mate_edit_distance_frac = 0.09f;

		//	max_no_mappings		= 32;		// Maximal no. of stored mappings
		//	mapping_counter_size = 1;
		max_no_mappings = 1 << 10;		// Maximal no. of stored mappings
		mapping_counter_size = 4;
		paired_reads = false;	// single-end reads

		sensitive_mode = true;
		sensitivity_factor = 2.5;	
	
		min_approx_indel_len = 3;		
		max_approx_indel_len = 50;		
		max_approx_indel_mismatches = (uint32_t) (max_no_errors * sensitivity_factor);

		high_confidence_sigmas = 4.0;
		penalty_saturation_sigmas = 7.0;
		clipping_distance_sigmas = 14.0;
		score_discretization_threshold = 0.5;
		hit_merging_threshold = 12;
		hit_merging_wrt_first = true;
		
		match_score = 1.0;
		mismatch_score = -5.0;
//		gap_open = -6.0; 
//		gap_extend = -1.0;
		clipping_score = -6;

		gap_ins_open = -5.0;
		gap_ins_extend = -0.4;
		gap_del_open = -5.0;
		gap_del_extend = -0.4;

		enable_boundary_clipping = true;
		enable_paired_clipping = true;
		min_clipped_factor = 1.0;

		mapq_mult = 90.0;
		mapq_div = 0.5;

		mask_low_quality_bases = 0;

		gzipped_SAM_level = 0;		// 0 - no compression
		store_BAM = false;
		use_stdout = false;

		verbosity_level = 1;
		filter = "";

		is_fasta = false;
		min_read_len = 135;
		max_read_len = 151;

		enable_mapping_indels = true;
		enable_var_indel_long = false;
		enable_var_indel_short = false;
		enable_var_snp = false;
		enable_short_reads = false;
		enable_short_indel_refinement = true;

		project_name = "whisper";	
		temp_prefix = "./whisper_temp_";
		read_group_line = "";
#ifdef _DEVELOPMENT_MODE
//		keep_temporary_files = true;
		keep_temporary_files = false;
#else
		keep_temporary_files = false;
#endif
	}
};

// ************************************************************************************
struct CParams : public CCmdParams
{
	bool ssd_mode;
	bool is_fasta;

	uint32_t no_thr_fastq_reader;
	uint32_t no_thr_splitter;
	uint32_t no_thr_mapping_cores;
	uint32_t no_thr_sam_generators;

	uint64_t max_mapping_memory;

	uint64_t max_res_dev_memory;
	uint64_t sam_buffer_memory;
	uint64_t sam_part_size;

	uint64_t block_overhead_size;
	uint64_t block_size;
	uint64_t no_fastq_blocks_per_thread;
	uint32_t max_fastq_rec_length;
	uint32_t max_no_errors;

	uint64_t bin_size;
	uint64_t read_len;
	uint64_t no_parts;
	uint64_t min_no_free_parts;
	uint64_t max_bin_part_size;

	uint64_t res_group_size;
	uint32_t no_res_groups;
	uint32_t no_res_parts;
	uint64_t min_no_free_group_parts;
	uint64_t max_group_part_size;

	uint32_t id_bits_total;
	uint32_t id_bits_subgroup;
	uint32_t id_bits_local;

	uint32_t max_cigar_len;
	uint32_t sa_prefix_overhead;
	bool constant_read_len;
	bool max_reads_compression;
	bool io_calls_in_serial_mode;

	vector<string> prefix_file_names;
	vector<uint32_t> prefix_map;

	vector<vector<pair<uint32_t, uint32_t>>> stage_segments;

	string read_group_id;

	CParams() : CCmdParams()
	{
		ssd_mode = false;

		io_calls_in_serial_mode = true;

		// Bits for parts of read ids
		id_bits_total = 40;		// total number of bits per id
		id_bits_subgroup = 3;		// number of bits for subgroups
		id_bits_local = 21;		// number of bits for id within a subgroup

		sa_prefix_overhead = 1;			// no. of extension of bin prefix for SA queries
	
		//	constant_read_len     = true;	// all reads are of the same length (so, read length is not stored for each read)
		constant_read_len = false;	// read lengths are stored (so, reads can be of different lenths)
		max_reads_compression = true;	// use 3-symbol into byte compression of reads (instead of default 2-in-1)
		//	max_reads_compression = false;	// use 3-symbol into byte compression of reads (instead of default 2-in-1)

		read_len = 10240 + 1;

		read_group_id = "";
	}

	CParams &operator=(const CCmdParams &cmd_params)
	{
		developer_mode = cmd_params.developer_mode;
		log_stats = cmd_params.log_stats;

		command_line = cmd_params.command_line;
		input_file_names = cmd_params.input_file_names;
		index_name = cmd_params.index_name;
		project_name = cmd_params.project_name;
		temp_prefix = cmd_params.temp_prefix;
		is_fasta = cmd_params.is_fasta;
		keep_temporary_files = cmd_params.keep_temporary_files;

		no_bins = cmd_params.no_bins;
		no_threads = cmd_params.no_threads;
		max_total_memory = cmd_params.max_total_memory;
		sa_in_ram = cmd_params.sa_in_ram;

		// Mapping 
		min_read_len = cmd_params.min_read_len;
		max_read_len = cmd_params.max_read_len;
		paired_reads = cmd_params.paired_reads;
		max_no_errors = cmd_params.max_no_errors;
		max_no_mappings = cmd_params.max_no_mappings;
		mapping_mode = cmd_params.mapping_mode;
		mapping_orientation = cmd_params.mapping_orientation;

		mapping_counter_size = cmd_params.mapping_counter_size;
		max_mate_distance = cmd_params.max_mate_distance;
		max_mate_edit_distance_frac = cmd_params.max_mate_edit_distance_frac;
		sensitive_mode = cmd_params.sensitive_mode;
		sensitivity_factor = cmd_params.sensitivity_factor;
		penalty_saturation_sigmas = cmd_params.penalty_saturation_sigmas;
		clipping_distance_sigmas = cmd_params.clipping_distance_sigmas;
		high_confidence_sigmas = cmd_params.high_confidence_sigmas;
		mask_low_quality_bases = cmd_params.mask_low_quality_bases;
		hit_merging_threshold = cmd_params.hit_merging_threshold;
		hit_merging_wrt_first = cmd_params.hit_merging_wrt_first;

		enable_boundary_clipping = cmd_params.enable_boundary_clipping;
		enable_paired_clipping = cmd_params.enable_paired_clipping;
		score_discretization_threshold = cmd_params.score_discretization_threshold;
		enable_short_indel_refinement = cmd_params.enable_short_indel_refinement;

		match_score = cmd_params.match_score;
		mismatch_score = cmd_params.mismatch_score;
//		gap_open = cmd_params.gap_open;
//		gap_extend = cmd_params.gap_extend;
		clipping_score = cmd_params.clipping_score;
		gap_ins_open = cmd_params.gap_ins_open;
		gap_ins_extend = cmd_params.gap_ins_extend;
		gap_del_open = cmd_params.gap_del_open;
		gap_del_extend = cmd_params.gap_del_extend;

		min_approx_indel_len = cmd_params.min_approx_indel_len;
		max_approx_indel_len = cmd_params.max_approx_indel_len;
		max_approx_indel_mismatches = cmd_params.max_approx_indel_mismatches;
		min_clipped_factor = cmd_params.min_clipped_factor;

		mapq_mult = cmd_params.mapq_mult;
		mapq_div = cmd_params.mapq_div;

		enable_mapping_indels = cmd_params.enable_mapping_indels;
		enable_var_indel_long = cmd_params.enable_var_indel_long;
		enable_var_indel_short = cmd_params.enable_var_indel_short;
		enable_var_snp = cmd_params.enable_var_snp;
		enable_short_reads = cmd_params.enable_short_reads;

		// Output format
		gzipped_SAM_level = cmd_params.gzipped_SAM_level;
		store_BAM = cmd_params.store_BAM;
		use_stdout = cmd_params.use_stdout;

		// User interface
		verbosity_level = cmd_params.verbosity_level;
		filter = cmd_params.filter;

		// parse RG line
		read_group_line = cmd_params.read_group_line;
		if (cmd_params.read_group_line != "") {
			// replace \\t with \t
			size_t index = 0;
			while ((index = read_group_line.find("\\t", index)) != string::npos) {
				read_group_line.replace(index, 2, "\t");
				index += 2;
			}

			std::istringstream iss(read_group_line);

			// extract id
			string token;
			while (iss >> token) {
				index = token.find("ID:", 0);
				if (index != string::npos) {
					read_group_id = token.substr(3);
				}
			}
		}

		return *this;
	}

	friend ostream& operator<<(ostream &out, const CParams &p);
};

// ************************************************************************************
struct CObjects
{
	CMemoryPool<uchar_t> *mp_fastq_blocks;
	CMemoryPool<uchar_t> *mp_reads;
	CMemoryPool<uchar_t> *mp_bins_read;
	CMemoryPool<uchar_t> *mp_bins_write;
	CMemoryPool<uchar_t> *mp_sam_parts;
	CMemoryPool<uchar_t> *mp_fastq_records;
	CMemoryPool<uchar_t> *mp_ext_cigar;
	CMemoryPool<uchar_t> *mp_cigar;
	CMemoryPool<uint32_t> *mp_cigar_bin;
	CMemoryPool<uchar_t> *mp_mdz;
	
	CRegisteringQueue<file_name_no_t> *q_file_names;
	CRegisteringQueue<fastq_block_t> *q_blocks;
	CRegisteringQueue<reads_bin_t> *q_bins_read;
	CRegisteringQueue<reads_bin_t> *q_bins_write;

	CRegisteringQueue<res_group_t> *q_map_res;
	CMemoryPool<uchar_t> *mp_map_res;
	
	// Ids of mapping result groups ready to convert to SAM
	CRegisteringQueue<uint32_t> *q_res_ids;

	// Mapped reads and raw reads for SAM generation
	CRegisteringQueue<mapped_reads_t> *q_map_reads;

	// SAM blocks
	CRegisteringQueue<sam_block_t> *q_sam_blocks;

	// Memory monitor
	CMemoryMonitor *mem_monitor;

	// Serial processing
	CSerialProcessing *serial_processing;

	// Pool of pointers
	CPtrPool *ptr_pool;

	// Reference
	CReference *reference;
	
	// Suffix arrays
	CSuffixArray *sa_dir;
	CSuffixArray *sa_rc;

#ifdef ENABLE_VCF_VARIANTS
	// Variants
	CVariantDB *variant_db;
#endif

	CRunningStats *running_stats;

	CProgress *progress;

	int64_t mem_alloc;
	mutex mtx;

	CObjects() : 
		mp_fastq_blocks(nullptr),
		mp_reads(nullptr),
		mp_bins_read(nullptr),
		mp_bins_write(nullptr),
		mp_sam_parts(nullptr),
		mp_fastq_records(nullptr),
		mp_ext_cigar(nullptr),
		mp_cigar(nullptr),
		mp_cigar_bin(nullptr),
		mp_mdz(nullptr),
		q_file_names(nullptr),
		q_blocks(nullptr),
		q_bins_read(nullptr),
		q_bins_write(nullptr),
		q_map_res(nullptr),
		mp_map_res(nullptr),
		q_res_ids(nullptr),
		q_map_reads(nullptr),
		q_sam_blocks(nullptr),
		mem_monitor(nullptr),
		serial_processing(nullptr),
		ptr_pool(nullptr),
		reference(nullptr),
		sa_dir(nullptr),
		sa_rc(nullptr),
#ifdef ENABLE_VCF_VARIANTS
		variant_db(nullptr),
#endif
		running_stats(nullptr),
		progress(nullptr),
		mem_alloc(0)
	{}

	void IncreaseMem(string s, int64_t x)
	{
		unique_lock<mutex> lck(mtx);
		mem_alloc += x;
	}

	void DecreaseMem(string s, int64_t x)
	{
		unique_lock<mutex> lck(mtx);
		mem_alloc -= x;
	}
};

#endif

// EOF
