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


#include "../common/joiner_mgr.h"

// ************************************************************************************
CJoinerMgr::CJoinerMgr(CParams *params, CObjects *objects)
{
	id_bits_total    = params->id_bits_total;
	id_bits_subgroup = params->id_bits_subgroup;
	id_bits_local    = params->id_bits_local;

	no_thr_fastq_reader = params->no_thr_fastq_reader;

	verbosity_level  = params->verbosity_level;

	q_map_reads      = objects->q_map_reads;
	q_res_ids        = objects->q_res_ids;

	completed = false;
}

// ************************************************************************************
CJoinerMgr::~CJoinerMgr()
{
}

// ************************************************************************************
// Insert FASTQ block into vector of blocks of a mapping bin
bool CJoinerMgr::PutFastqBlock(fastq_block_t fastq_block, bool single_end)
{
	lock_guard<mutex> lck(mtx);

	uint32_t bin_id = (uint32_t) (fastq_block.id_range >> (id_bits_local + id_bits_subgroup)); 

	auto p = data.find(bin_id);
	if(p == data.end())
	{
		p = data.insert(make_pair(bin_id, mapped_reads_t())).first;
		p->second.single_end = single_end;
	}

	p->second.fastq_blocks.push_back(fastq_block);

	return true;
}

// ************************************************************************************
// Insert mapping group results 
bool CJoinerMgr::PutResultsBin(uint32_t bin_id, results_bin_t results_bin)
{
	lock_guard<mutex> lck(mtx);

	auto p = data.find(bin_id);
	if(p == data.end())
		p = data.insert(make_pair(bin_id, mapped_reads_t())).first;

	p->second.results = results_bin;

	// If all FASTQ blocks and results ready, then send for SAM generation
	if(p->second.all_blocks)
		send_bin_for_processing(bin_id);
	
	return true;
}

// ************************************************************************************
// Mark that all FASTQ blocks for a results bin are inserted
void CJoinerMgr::MarkAllBinBlocks(uint32_t bin_id)
{
	lock_guard<mutex> lck(mtx);

	auto p = data.find(bin_id);

	p->second.all_blocks = true;

	// If all FASTQ blocks and results ready, then send for SAM generation
	if(p->second.results.size)
		send_bin_for_processing(bin_id);
}

// ************************************************************************************
void CJoinerMgr::MarkBinsCompleted()
{
	lock_guard<mutex> lck(mtx);

	if(no_thr_fastq_reader == 0)
		complete();
	else
		completed = true;
}

// ************************************************************************************
// Mark that FASTQ reader finished reading 
void CJoinerMgr::MarkFastqReaderCompleted()
{
	lock_guard<mutex> lck(mtx);

	if(--no_thr_fastq_reader == 0 && completed)
		complete();
}

// ************************************************************************************
// Send bin for processing - insert bin into a queue 
void CJoinerMgr::send_bin_for_processing(uint32_t bin_id)
{
	if(verbosity_level >= 3)
		cerr << "JoinerMgr: sending result group for converting into SAM: " << bin_id << " size of data: " << data.size() << "\n";

	q_map_reads->Push(data[bin_id]);

	data.erase(bin_id);

	if(data.empty() && completed && no_thr_fastq_reader == 0)
		complete();
}

// ************************************************************************************
// Marks that JoinerMgr finished processing the results
void CJoinerMgr::complete()
{
	if(verbosity_level >= 3)
		cerr << "JoinerMgr: mapping results reading finished\n";

	q_map_reads->MarkCompleted();
}

// EOF
