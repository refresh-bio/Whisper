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


#ifndef _MAPPER_H
#define _MAPPER_H

#include <vector>
#include <string>
#include <thread>
#include <algorithm>
#include <utility>

#include "../common/idstore.h"
#include "../common/utils.h"
#include "../common/mmgr.h"
#include "../common/queue.h"
#include "../common/fastq_reader.h"
#include "../common/stats.h"
#include "../common/sa.h"
#include "../common/joiner_mgr.h"
#include "ref_desc.h"
#include "reads.h"
#include "bins.h"
#include "mapping_core.h"
#include "../common/params.h"
#include "sam.h"
#include "sam_writer.h"

using namespace std;

class CMapper
{
	CParams params;
	CObjects objects;
	CRefSeqDesc ref_seq_desc;
	
	int64_t total_input_file_size;

	bool adjust_threads(bool single_thr_readers);
	bool adjust_bins();
	bool adjust_memory_splitting();
	bool adjust_memory_mapping();
	bool adjust_memory_postprocessing();
	bool adjust_stage_segments();
	bool adjust_max_and_min_read_len(vector<uint64_t> &hist_read_len);
	
	bool prepare_running_stats(bool before_preprocessing);

	bool prepare_reference();
	bool prepare_reference_sa();
	bool release_reference();
	bool release_reference_sa();
	
	bool prepare_serial_processing(bool call_in_serial_mode);
	bool release_serial_processing();

	uint32_t get_next_stage(uint32_t c_stage, uint32_t m_stage);

	// Mapping stages
	bool reads_splitting();						// Initial reads splitting
	bool reads_mapping_first_stratum();			// Reads mapping (first stratum)
	bool reads_mapping_second_stratum();		// Reads mapping (second stratum)
	bool reads_mapping_all_strata();			// Reads mapping (all strata)
	bool stage(uint32_t stage_major, uint32_t stage_minor, uint32_t prev_stage_major, uint32_t prev_stage_minor, mapping_mode_t mapping_mode, bool sensitive_mode = false, bool final_substage = false);

	bool reads_mapping_single_stage(uint32_t stage_major, uint32_t stage_minor, mapping_mode_t mapping_mode, bool sensitive_mode);
	bool reads_mapping_stage_range(uint32_t stage_major_from, uint32_t stage_minor_from, uint32_t stage_major_to, uint32_t stage_minor_to, mapping_mode_t mapping_mode, bool sensitive_mode);
	bool reads_postprocessing(mapping_mode_t mapping_mode);

	bool check_input_files();
	bool load_ref_seq_desc();
	bool prepare_output_files(vector<string> &header_SAM, vector<pair<string, uint32_t>> &header_BAM);

public:
	CMapper();
	~CMapper();

	bool SetRunningStats(CRunningStats *_running_stats);
	bool SetParams(const CCmdParams &cmd_params);

	bool StartMapping();

	// Only in development mode
	bool StartMapping(uint32_t stage_major, uint32_t stage_minor, bool sensitive_mode);
	bool StartMapping(uint32_t stage_major_from, uint32_t stage_minor_from, uint32_t stage_major_to, uint32_t stage_minor_to, bool sensitive_mode);
	bool StartPostProcessing();
};

#endif

// EOF
