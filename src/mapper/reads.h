// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.0
// Date    : 2017-11-30
// License : GNU GPL 3
// *******************************************************************************************


#ifndef _READS_H
#define _READS_H

#include "../common/defs.h"
#include "../common/utils.h"
#include "../common/mmgr.h"
#include "../common/queue.h"
#include "../common/fastq_reader.h"
#include "../common/stats.h"
#include "../common/timer.h"
#include "../common/params.h"

// ************************************************************************************
//
// ************************************************************************************
class CReadsCollector
{
	CMemoryPool<uchar_t> *mp_bins_write;
	CRegisteringQueue<reads_bin_t> *q_bins_write; 
	uint32_t no_bins;
	vector<uint32_t> prefix_map;
	CRunningStats *running_stats;
	CStopWatch watch;

	uchar_t** buffers;
	uint32_t* buffer_sizes;
	uint32_t* read_counters;
	uint32_t part_size;
	uint32_t read_id_bytes;
	bool constant_read_len;
	bool max_reads_compression;

	void prepare_buffers();
	inline uint32_t pack_read(uchar_t* dest, read_id_t id, uchar_t *data, uint32_t size, uint32_t best_error);
	inline uint32_t copy_read(uchar_t* dest, read_id_t id, uchar_t *data, uint32_t size, uint32_t best_error);

public:
	CReadsCollector(CParams *params, CObjects *objects);
	~CReadsCollector();

	bool Push(uint32_t bin_id, read_id_t id, uchar_t *data, uint32_t size, uint32_t best_error, bool is_packed = true);
	bool Complete();
};

// ************************************************************************************
//
// ************************************************************************************
class CReadsDeliverer
{
	CRunningStats *running_stats;
	uint32_t stage_id;
	CStopWatch watch;

	CMemoryPool<uchar_t> *mp_reads;

	uint16_t decompress_lut[128];

	bool is_valid;
	bool is_sorted;
	bool is_compressed;

	uint32_t min_read_len;
	uint32_t max_read_len;
	bool constant_read_len;
	bool max_reads_compression;

	reads_bin_t bin;
	uint64_t pos;
	uint32_t read_id_len;
	uint32_t start_pos;			// position of the starting symbol for sorting purposes in read
	uint32_t end_pos;			// position of the ending symbol for sorting purposes

	uchar_t **sorted_ptrs;
	uchar_t *buf_data_sft;

	void sort_reads();
	void decompress_reads();
	void prepare_lut();

public:
	CReadsDeliverer(CParams *params, uint32_t _read_id_len, uint32_t _start_pos, uint32_t _end_pos, CRunningStats *_running_stats, 
		CMemoryPool<uchar_t> *_mp_reads, uint32_t _stage_id);
	~CReadsDeliverer();

	bool SetBin(reads_bin_t &_bin);
	uint64_t GetNoReads();
	bool Start();
	bool Pop(read_id_t &id, uchar_t* &data, uchar_t* &data_sft, uint32_t &size, uint32_t &best_error);
};


// ************************************************************************************
//
// ************************************************************************************
class CReadsSplitter
{
	CMemoryPool<uchar_t> *mp_bins_write;
	CMemoryPool<uchar_t> *mp_reads;
	CRegisteringQueue<reads_bin_t> *q_bins_write;
	uint32_t no_bins;
	uint32_t bin_prefix;
	uint32_t read_id_len;
	bool is_fasta;
	uint32_t mask_low_quality_bases;

	vector<uint32_t> prefix_map;
	uint32_t next_start_pos;
	uint32_t next_end_pos;

	uchar_t codes[256];
	uchar_t *block;
	uint32_t block_size, block_pos;

	uint64_t n_reads;
	vector<uint64_t> hist_read_len;

	CReadsCollector *reads_collector;
	CRunningStats *running_stats;
	CStopWatch watch;

	bool get_seq(uchar_t *seq, uint32_t &seq_size);

public:
	CReadsSplitter(CParams *params, CObjects *objects, bool _is_fasta, uint32_t first_major_stage);
	~CReadsSplitter();

	bool Process(uchar_t *_block, uint32_t _block_size, read_id_t id_range);
	void Complete();
	void GetHistReadLen(vector<uint64_t> &_hist_read_len)
	{
		_hist_read_len = hist_read_len;
	}
};

// ************************************************************************************
//
// ************************************************************************************
class CWReadsSplitter
{
	CReadsSplitter *rs;

	CRegisteringQueue<fastq_block_t> *q_fastq_blocks;
	CRegisteringQueue<reads_bin_t> *q_bins_write;
	CMemoryPool<uchar_t> *mp_blocks;

	vector<uint64_t> hist_read_len;
	uint32_t verbosity_level;

public:
	CWReadsSplitter(CParams *params, CObjects *objects, bool _is_fasta, uint32_t first_major_stage);
	~CWReadsSplitter();

	CWReadsSplitter(const CWReadsSplitter &x) = delete;
	CWReadsSplitter(CWReadsSplitter &&x) = delete;

	void operator()();

	void GetHistReadLen(vector<uint64_t> &_hist_read_len)
	{
		_hist_read_len = hist_read_len;
	}
};

// ************************************************************************************
//
// ************************************************************************************
class CReadsReader
{
	uint32_t read_id_len;
	bool is_fasta;

	uint64_t n_reads;

	uchar_t *block;
	uint64_t block_size;
	uint64_t block_pos;

	uchar_t *id;
	uchar_t *sequence;
	uchar_t *plus;
	uchar_t *quality;

	CRunningStats *running_stats;
	CStopWatch watch;

	CMemoryPool<uchar_t> *mp_fastq_records;

	bool get_seq(uchar_t *id, uchar_t *sequence, uchar_t *plus, uchar_t *quality);

public:
	CReadsReader(CParams *params, CObjects *objects, bool _is_fasta);
	~CReadsReader();

	bool SetBlock(uchar_t *_block, uint64_t _block_size);
	bool Pop(uchar_t* &_id, uchar_t* &_sequence, uchar_t* &_plus, uchar_t* &_quality,
		uint32_t &id_len, uint32_t &sequence_len, uint32_t &plus_len, uint32_t &quality_len);
};

#endif

// EOF
