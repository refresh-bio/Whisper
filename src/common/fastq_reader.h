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


#ifndef _FASTQ_READER_H
#define _FASTQ_READER_H

#include "../common/defs.h"
#include "../common/mmgr.h"
#include "../common/queue.h"
#include "../common/idstore.h"
#include "../common/params.h"
#include "../common/joiner_mgr.h"
#include "../mapper/simd_utils.h"
#include <cstdio>
#include <iostream>
#include <string>

#ifndef _DEBUG
#include "../libs/zlib.h"
#include "../libs/bzlib.h"
#endif

using namespace std;

enum class FASTQ_reader_thr_state_t {wait, can_read, data_ready, eof, complete};

// ************************************************************************************
// Wrapper for zlib
// ************************************************************************************
class CGZFile {
	CSerialProcessing *serial_processing;
	int buffer_size;
	uchar_t *buffer;
	FILE *in;
	size_t in_size;
	size_t in_buffer_size;
	bool is_eof;

#ifndef _DEBUG
	z_stream stream;
#endif

public:
	CGZFile(int _buffer_size, CSerialProcessing *_serial_processing = nullptr) : 
		serial_processing(_serial_processing), buffer_size(_buffer_size)
	{
		buffer = new uchar_t[buffer_size];
		in = nullptr;
	}

	~CGZFile()
	{
		if (in)
			Close();

		if (buffer)
			delete[] buffer;
	}

	bool Open(string file_name);
	bool Close();
	bool Read(uchar_t *ptr, uint32_t to_read, uint32_t &readed, uint32_t &raw_readed);
	bool Eof(void);
};

// ************************************************************************************
// FASTA/FASTQ reader class
// ************************************************************************************
class CFastqReader {
protected:
	typedef enum class mode_t_ {m_plain, m_gzip, m_bzip2} mode_t;

	CMemoryPool<uchar_t> *mem_pool;
	CProgress *progress;

	string input_file_name;
	bool is_fasta;
	mode_t mode;

	FILE *in;
#ifndef _DEBUG
	gzFile_s *in_gzip;
	BZFILE *in_bzip2;
	int bzerror;

	CGZFile *in_gz_file;
#endif

	uint32_t part_size;
	
	uchar_t *internal_buffer;
	uint32_t part_filled;
	uint32_t part_no_reads;
	uint32_t last_read_start_pos;
	uint32_t total_filled;

	CSerialProcessing *serial_processing;

#ifndef _DEBUG
	uint32_t gzip_buffer_size;
	uint32_t bzip2_buffer_size;
#endif

	size_t(*ptr_CountEOLs)(uchar_t *, size_t);

	bool SkipNextEOL(uchar_t *part, uint32_t &pos, uint32_t max_pos);
	bool Eof();
	void close_files();
	void count_reads(uchar_t *part, uint32_t filled, uint32_t &no_reads, uint32_t &last_read_start_pos);

	uint32_t find_nth_read(uchar_t *part, uint32_t size, uint32_t no_reads);

public:
	CFastqReader(bool _is_fasta, uint32_t _part_size, CObjects *_objects);
	~CFastqReader();

	static uint32_t OVERHEAD_SIZE;

	void Restart();
	bool Open(string _input_file_name);
	bool GetPartInfo(uchar_t *&part, uint32_t &size, uint32_t &no_reads);
	bool GetPart(uchar_t *part, uint32_t &size, uint32_t no_reads);
	void DisableProgress(void)
	{
		progress = nullptr;
	}
};

// ************************************************************************************
// FASTA/FASTQ reader thread class for pre/post-processing stage (base class)
// ************************************************************************************
class CPrePostBaseFastqReader {
protected:
	CFastqReader *fqr1, *fqr2;

	CMemoryPool<uchar_t> *mem_pool;
	CRegisteringQueue<file_name_no_t> *q_file_names;
	CIDStore *id_store;
	CRunningStats *running_stats;

	uint32_t verbosity_level;

	void load_single_file(file_name_no_t file_name_no);
	void load_pair_files(file_name_no_t file_name_no);
	void load_pair_files_thr(file_name_no_t file_name_no);

	virtual void process_single_block(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id, uchar_t* part, uint32_t part_filled) = 0;
	virtual void process_pair_blocks(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id1, uchar_t* part1, uint32_t part_filled1, 
		read_id_t &id2, uchar_t* part2, uint32_t part_filled2) = 0;
	virtual void store_stats(double time) = 0;
	virtual void mark_eof(file_name_no_t file_name_no, read_id_t id) = 0;
	virtual void end_of_thread() = 0;

public:
	CPrePostBaseFastqReader(bool _is_fasta, CParams *params, CObjects *objects, CIDStore *_id_store);
	virtual ~CPrePostBaseFastqReader();

	void operator()();
};

// ************************************************************************************
// FASTA/FASTQ reader thread class for preprocessing stage
// ************************************************************************************
class CPreFastqReader : public CPrePostBaseFastqReader {
	CRegisteringQueue<fastq_block_t> *q_blocks;

	void process_single_block(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id, uchar_t* part, uint32_t part_filled);
	void process_pair_blocks(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id1, uchar_t* part1, uint32_t part_filled1, 
		read_id_t &id2, uchar_t* part2, uint32_t part_filled2);
	void store_stats(double time);
	void mark_eof(file_name_no_t file_name_no, read_id_t id) {};
	void end_of_thread();

public:
	CPreFastqReader(bool _is_fasta, CParams *params, CObjects *objects, CIDStore *_id_store) :
		CPrePostBaseFastqReader(_is_fasta, params, objects, _id_store)
	{
		q_blocks     = objects->q_blocks;
	};
	~CPreFastqReader() {};
};

// ************************************************************************************
// FASTA/FASTQ reader thread class for postprocessing stage
// ************************************************************************************
class CPostFastqReader : public CPrePostBaseFastqReader {
	CRegisteringQueue<uint32_t> *q_res_ids;
	CJoinerMgr *joiner_mgr;
	uint32_t id_bits_local;
	uint32_t id_bits_subgroup;
	read_id_t sub_block_mask;

	int64_t send_bytes;

	void process_single_block(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id, uchar_t* part, uint32_t part_filled);
	void process_pair_blocks(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id1, uchar_t* part1, uint32_t part_filled1, 
		read_id_t &id2, uchar_t* part2, uint32_t part_filled2);
	void store_stats(double time);
	void mark_eof(file_name_no_t file_name_no, read_id_t id);
	bool is_new_group(read_id_t prev_id, read_id_t new_id);
	void end_of_thread();

public:
	CPostFastqReader(bool _is_fasta, CParams *params, CObjects *objects, CIDStore *_id_store, CJoinerMgr *_joiner_mgr) :
		CPrePostBaseFastqReader(_is_fasta, params, objects, _id_store)
	{
		q_res_ids        = objects->q_res_ids;
		joiner_mgr       = _joiner_mgr;
		id_bits_local    = params->id_bits_local;
		id_bits_subgroup = params->id_bits_subgroup;
		sub_block_mask   = (((read_id_t) 1) << id_bits_subgroup) - 1;

		send_bytes = 0;

		if (fqr1)
			fqr1->DisableProgress();
		if (fqr2)
			fqr2->DisableProgress();
	};
	~CPostFastqReader() 
	{
		running_stats->AddValues(1000001, send_bytes);
	};
};

#endif

// EOF
