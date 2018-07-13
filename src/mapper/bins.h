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


#ifndef _BINS_H
#define _BINS_H

#include "../common/defs.h"
#include "../common/utils.h"
#include "../common/mmgr.h"
#include "../common/queue.h"
#include "../common/stats.h"
#include "../common/timer.h"
#include "../common/params.h"
#include "reads.h"
#include <memory>
#include <list>
#include <vector>

// ************************************************************************************
//
class CBinsWriter
{
	uint32_t no_bins;
	CMemoryPool<uchar_t> *mp_bins;
	CRegisteringQueue<reads_bin_t> *q_bins_write;
	vector<shared_ptr<CMapperFile>> f_bins;
	CRunningStats *running_stats;
	uint32_t stage_id;
	CStopWatch watch;
	uint32_t verbosity_level;
	bool keep_temporary_files;
	bool final_substage;

	string directory;
	string name_prefix;
	uint64_t min_no_free_parts;
	uint64_t max_mem_single_file;
	uint64_t largest_bin_size;
	uint32_t largest_bin_id;
	uint32_t n_parts;

	vector<string> prefix_file_names;

	vector<uint64_t> buffer_sizes;
	vector<uint64_t> read_counters;
	vector<vector<reads_bin_t>> buffer_parts;
	uchar_t *writing_buffer;

	uint64_t n_total_parts;
	uint64_t n_write_bytes;

	string get_name(int n);
	bool open_files();
	bool close_files();
	bool write_bin(uint32_t bin_id);
	void store_stats();

public:
	CBinsWriter(CParams *params, CObjects *objects, string _name_prefix, uint32_t _stage_id, bool _final_substage);
	~CBinsWriter();

	void operator()();
};


// ************************************************************************************
//
class CBinsReader
{
	string directory;
	string name_prefix;
	vector<string> prefix_file_names;
	CMemoryMonitor *mem_monitor;
	CRegisteringQueue<reads_bin_t> *q_bins_read;
	vector<uint64_t> read_counters;
	uint32_t no_bins;
	CRunningStats *running_stats;
	CStopWatch watch;
	uint32_t stage_id;
	uint32_t verbosity_level;
	bool max_reads_compression;
	uint32_t read_len;
	bool keep_temporary_files;

	string get_name(int n);
	bool is_mono_symbol(int n);
	bool load_stats();

public:
	CBinsReader(CParams *params, CObjects *objects, string _name_prefix, uint32_t _stage_id);
	~CBinsReader();

	void operator()();
};


// ************************************************************************************
//
class CBinPrefixes
{
	bool is_valid;

	uint32_t max_prefix_len;
	uint32_t prefix_table_size;
	uint32_t no_bins;
	uint32_t verbosity_level;

	char codes[4];
	uint32_t decodes[256];

	vector<uint32_t> prefix_map;			// maps long prefixes onto bin ids
	vector<string> prefix_file_names;		// file names of prefixes
	vector<uint64_t> small_lut;				// counts of prefix (of length max_prefix_len) appearance in SA

	inline uint64_t count_suffixes_from(string str);

public:
	CBinPrefixes(uint32_t _max_prefix_len, uint32_t _verbosity_level);
	~CBinPrefixes();

	void SetLUT(uint32_t *lut_data, uint32_t lut_prefix_len);
	bool ReduceMap(uint32_t _no_bins);
	bool GetMaps(vector<uint32_t> &_prefix_map, vector<string> &_prefix_file_names);
};

#endif

// EOF
