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


#ifndef _SAM_WRITER_H
#define _SAM_WRITER_H

#include "../common/defs.h"
#include "../common/queue.h"
#include "../common/mmgr.h"
#include "../common/stats.h"
#include "../common/types.h"
#include "../common/params.h"
#include "reads.h"
#include "../common/joiner_mgr.h"
#include "soft_clipping.h"
#include "ref_desc.h"
#include "LevMyers.h"

// ************************************************************************************
class CSamWriter
{
	CMapperFile file_mapped;
	vector<sam_block_t> blocks;
	uint64_t total_size_mapped;
	uint64_t max_size;
	uint64_t writes_no;
	uint32_t verbosity_level;
	uint32_t gzipped_SAM_level;
	bool store_BAM;
	bool use_stdout;

	CRegisteringQueue<sam_block_t> *q_sam_blocks;

	CMemoryPool<uchar_t> *mp_sam_parts;
	uchar_t *buffer_mapped;	

	CSerialProcessing *serial_processing;

	string gz_raw_buffer;
	CGzipMember *gzip;

	mutable mutex mtx;								// The mutex to synchronise on

	void write(CMapperFile &file, uint64_t &total_size, uchar_t* buffer, string s);
	void gz_write(CMapperFile &file, uint64_t &total_size, uchar_t* buffer, string s);
	void gz_begin_member(CMapperFile &file, uint64_t &total_size, uchar_t *buffer);
	void gz_end_member(CMapperFile &file, uint64_t &total_size, uchar_t *buffer);
	void complete();

	void write_SAM_header(vector<string> &header_SAM);
	void write_BAM_header(vector<string> &header_SAM, vector<pair<string, uint32_t>> &header_BAM);

public:
	CSamWriter(CParams *params, CObjects *objects, vector<string> &header_SAM, vector<pair<string, uint32_t>> &header_BAM);
	~CSamWriter();

	void operator()();
};

#endif

// EOF
