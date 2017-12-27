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


#ifndef _RESULTS_H
#define _RESULTS_H

#include "../common/defs.h"
#include "../common/idstore.h"
#include "../common/queue.h"
#include "../common/mmgr.h"
#include "../common/stats.h"
#include "../common/types.h"
#include "../common/params.h"
#include "../common/joiner_mgr.h"

#include <vector>
#include <algorithm>

using namespace std;


// ************************************************************************************
class CMappingResultsCollector
{
	CRunningStats *running_stats;
	uint32_t stage_id;

	uint32_t no_res_groups;
	CMemoryPool<uchar_t> *mp_map_res;
	CRegisteringQueue<res_group_t> *q_map_res;
	uint32_t res_group_size;

	uint32_t total_bits;
	uint32_t group_bits;
	uint32_t id_bytes;
	uint32_t max_no_mappings;
	uint32_t mapping_counter_size;
	uint32_t max_counter_value;
	uint32_t buffer_reserve;

	uchar_t **buffers;
	uint32_t *buffer_sizes;

	read_id_t prev_read_id;
	uint32_t cur_read_count;
	uint32_t true_read_count;

	uint64_t push_counter;
	uint64_t push_unique_counter;
	uint64_t send_bytes;
	uint64_t added_bytes;

public:
	CMappingResultsCollector(CParams *params, CObjects *objects, uint32_t _stage_id);
	~CMappingResultsCollector();

	void Push(read_id_t id, ref_pos_t pos, genome_t direction, uint32_t no_differences);
	void Complete();
};

// ************************************************************************************
class CMappingResultsDeliverer
{
	CRunningStats *running_stats;

	string directory;
	string prefix_name;
	shared_ptr<CMapperFile> f_group;

	CRegisteringQueue<uint32_t> *q_res_ids;
	CJoinerMgr *joiner_mgr;
	CMemoryMonitor *mem_monitor;
	CSerialProcessing *serial_processing;

	bool keep_temporary_files;

	uint32_t verbosity_level;
	
	string get_name(int n);
	
public:
	CMappingResultsDeliverer(CParams *params, CObjects *objects, string _prefix_name, CJoinerMgr *_joiner_mgr);
	~CMappingResultsDeliverer();

	void operator()();
};

// ************************************************************************************
class CResultGroupsWriter
{
	CRunningStats *running_stats;

	string directory;
	string prefix_name;
	vector<shared_ptr<CMapperFile>> f_groups;

	CMemoryPool<uchar_t> *mp_map_res;
	CRegisteringQueue<res_group_t> *q_map_res;
	uint32_t no_groups;
	uint64_t min_no_free_parts;
	uint64_t max_mem_single_file;

	uint64_t largest_group_size;
	uint32_t largest_group_id;

	vector<uint64_t> buffer_sizes;
	vector<vector<res_group_t>> buffer_parts;
	uchar_t *writing_buffer;

	uint64_t n_total_parts;
	uint64_t n_write_bytes;
	uint64_t n_parts;

	string get_name(int n);
	
	bool open_files();
	bool close_files();
	bool write_part(uint32_t group_id);

public:
	CResultGroupsWriter(CParams *params, CObjects *objects, string _prefix_name);
	~CResultGroupsWriter();

	void operator()();
};

#endif

// EOF
