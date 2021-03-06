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


#ifndef _JOINER_MGR_H
#define _JOINER_MGR_H

#include <thread>

#include "../common/defs.h"
#include "idstore.h"
#include "queue.h"
#include "params.h"
#include "../common/types.h"

#include <map>

using namespace std;

class CJoinerMgr
{
	uint32_t id_bits_total;
	uint32_t id_bits_subgroup;
	uint32_t id_bits_local;

	uint32_t no_thr_fastq_reader;

	map<uint32_t, mapped_reads_t> data;
	bool completed;
	uint32_t verbosity_level;
	CRegisteringQueue<mapped_reads_t> *q_map_reads;
	CRegisteringQueue<uint32_t> *q_res_ids;

	mutable mutex mtx;								

	void send_bin_for_processing(uint32_t bin_id);
	void complete();

public:
	CJoinerMgr(CParams *params, CObjects *objects);
	~CJoinerMgr();

	bool PutFastqBlock(fastq_block_t fastq_block, bool single_end);
	bool PutResultsBin(uint32_t bin_id, results_bin_t results_bin);
	void MarkAllBinBlocks(uint32_t bin_id);

	void MarkFastqReaderCompleted();
	void MarkBinsCompleted();
};

#endif

// EOF
