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


#include "results.h"
#include <algorithm>
#include <utility>
#include "../libs/asmlib.h"
#include "../common/timer.h"
#include "../common/defs.h"


// ************************************************************************************
// CMappingsHeapGatherer
// ************************************************************************************

// ************************************************************************************
CMappingsHeapGatherer::CMappingsHeapGatherer(size_t _max_size)
{
	Clear(_max_size);
}

// ************************************************************************************
void CMappingsHeapGatherer::Clear(size_t _max_size)
{
	max_size = _max_size;
	v_mappings.resize(max_size);
	size = 0;
}

// ************************************************************************************
void CMappingsHeapGatherer::Push(uint32_t pos, genome_t direction, uint32_t no_errors)
{
	uint64_t x;

	x = ((uint64_t)no_errors) << 40;
	x += ((uint64_t)pos) << 8;
	x += (uint64_t) direction;
	
	if (size < max_size)
	{
		v_mappings[size++] = x;
		push_heap(v_mappings.begin(), v_mappings.begin() + size);
	}
	else
	{
		if (x < v_mappings.front())
		{
			pop_heap(v_mappings.begin(), v_mappings.end());
			v_mappings.back() = x;
			push_heap(v_mappings.begin(), v_mappings.end());
		}
	}
}

// ************************************************************************************
bool CMappingsHeapGatherer::Empty()
{
	return size == 0;
}

// ************************************************************************************
uint64_t CMappingsHeapGatherer::Pop()
{
	uint64_t r = 0;
		
	if (size)
	{
		r = v_mappings.front();
		pop_heap(v_mappings.begin(), v_mappings.begin() + size);
		--size;
	}

	return r;
}

// ************************************************************************************
uint64_t CMappingsHeapGatherer::PopUnsorted()
{
	return v_mappings[--size];
}

// ************************************************************************************
uint32_t CMappingsHeapGatherer::DecodePos(uint64_t x) const
{
	return (uint32_t)((x >> 8) & 0xffffffffull);
}

// ************************************************************************************
uint32_t CMappingsHeapGatherer::DecodeNoErrors(uint64_t x) const
{
	return (uint32_t)(x >> 40);
}

// ************************************************************************************
genome_t CMappingsHeapGatherer::DecodeDir(uint64_t x) const
{
	return (genome_t)(x & 0xffull);
}


// ************************************************************************************
// CMappingResultsCollector
// ************************************************************************************

// ************************************************************************************
CMappingResultsCollector::CMappingResultsCollector(CParams *params, CObjects *objects, uint32_t _stage_id)
{
	running_stats = objects->running_stats;
	stage_id = _stage_id;

	q_map_res      = objects->q_map_res;
	mp_map_res     = objects->mp_map_res;
	no_res_groups  = params->no_res_groups;
	res_group_size = (uint32_t) params->res_group_size;

	total_bits = params->id_bits_total;
	group_bits = params->id_bits_subgroup + params->id_bits_local;

	max_no_mappings      = params->max_no_mappings;
	mapping_counter_size = params->mapping_counter_size;  
	max_counter_value    = (uint32_t) ((1ull << (mapping_counter_size*8)) - 1);

	id_bytes   = BITS2BYTES(total_bits);

	buffer_reserve = id_bytes + mapping_counter_size + sizeof(ref_pos_t) + 2 + 1;

	buffers      = new uchar_t*[no_res_groups];
	buffer_sizes = new uint32_t[no_res_groups];

	for(uint32_t i = 0; i < no_res_groups; ++i)
	{
		mp_map_res->Reserve(buffers[i]);
		buffer_sizes[i] = 0;
	}

	prev_read_id = empty_read_id;

	push_counter = 0;
	push_unique_counter = 0;
	sent_bytes = 0;
	added_bytes = 0;
}

// ************************************************************************************
CMappingResultsCollector::~CMappingResultsCollector()
{
	if(buffers)
	{
		for(uint32_t i = 0; i < no_res_groups; ++i)
		{
			if(buffer_sizes[i])
			{
				cerr << "!!! buffer_sizes: " << buffer_sizes[i] << "\n";
			}
			if(buffers[i])
				mp_map_res->Free(buffers[i]);
		}

		delete[] buffers;
	}

	if(buffer_sizes)
		delete[] buffer_sizes;
}

// ************************************************************************************
// Store description of the mapping of a read:
// * read_id		- (default: 5B)
// * count			- (fixed: 1B)
// * data: (fixed: 6B x count, but not more than some threshold)
//   - pos_in_ref	- (fixed: 4B)
//   - direction	- (fixed: 1B)
//   - errors		- (fixed: 1B)
void CMappingResultsCollector::Push(read_id_t id, ref_pos_t pos, genome_t direction, uint32_t no_differences)
{
	uint32_t group_id = (uint32_t) (id >> group_bits);

	if(group_id >= no_res_groups)
	{
		cerr << "Error: " << group_id << "  " << no_res_groups << "\n";
		exit(1);
	}

	++push_counter;

	// If the group buffer is almost full, it is sent to queue for storage
	if(buffer_sizes[group_id] + buffer_reserve >= res_group_size)
	{
		q_map_res->Push(res_group_t(group_id, buffers[group_id], buffer_sizes[group_id]));
		sent_bytes += buffer_sizes[group_id];
		mp_map_res->Reserve(buffers[group_id]);
		buffer_sizes[group_id] = 0;
		
		prev_read_id    = empty_read_id;
		cur_read_count  = 0;
		true_read_count = 0;
	}

	if(id != prev_read_id)			// new read
	{
		StoreUInt(buffers[group_id]+buffer_sizes[group_id], id, id_bytes);			// id
		buffer_sizes[group_id] += id_bytes;
		added_bytes += id_bytes;
//		buffers[group_id][buffer_sizes[group_id]++] = 0;							// counter
		StoreUInt(buffers[group_id]+buffer_sizes[group_id], 0, mapping_counter_size);
		buffer_sizes[group_id] += mapping_counter_size;
		added_bytes += mapping_counter_size;
		prev_read_id    = id;
		cur_read_count  = 0;
		true_read_count = 0;

		push_unique_counter++;
	}

	uint32_t rec_size = sizeof(ref_pos_t) + 1 + 1;

	// !!! Tymczasowo ignorujemy max_no_mappings
/*	if(true_read_count >= max_no_mappings)		// do not store explicitly the mappings for reads that maps in too many places
	{
		if(true_read_count < max_counter_value)
			IncrementUInt(buffers[group_id]+(buffer_sizes[group_id] - mapping_counter_size - max_no_mappings*rec_size), mapping_counter_size);
//			++buffers[group_id][buffer_sizes[group_id]-1 - max_no_mappings*rec_size];		// increment counter
		++true_read_count;

		return;
	}
	*/

//	++buffers[group_id][buffer_sizes[group_id]-1 - cur_read_count*rec_size];		// increment counter
	IncrementUInt(buffers[group_id]+(buffer_sizes[group_id] - mapping_counter_size - cur_read_count*rec_size), mapping_counter_size);

	StoreUInt(buffers[group_id]+buffer_sizes[group_id], pos, sizeof(ref_pos_t));
	buffer_sizes[group_id] += sizeof(ref_pos_t);
	buffers[group_id][buffer_sizes[group_id]++] = (uchar_t) direction;
	buffers[group_id][buffer_sizes[group_id]++] = no_differences;
	added_bytes += sizeof(ref_pos_t) + 2;
	++cur_read_count;
	++true_read_count;
}

// ************************************************************************************
void CMappingResultsCollector::Complete()
{
	for(uint32_t i = 0; i < no_res_groups; ++i)
	{
		if(buffer_sizes[i])
		{
			q_map_res->Push(res_group_t(i, buffers[i], buffer_sizes[i]));
			sent_bytes += buffer_sizes[i];
			mp_map_res->Reserve(buffers[i]);
		}
		buffer_sizes[i] = 0;
	}

	running_stats->AddValues(STAT_PUSH_RESULTS    , (int64_t) push_counter);
	running_stats->AddValues(STAT_PUSH_RESULTS_UNQ, (int64_t) push_unique_counter);
	running_stats->AddValues(STAT_PUSH_RESULTS_SEND_BYTES, (int64_t) sent_bytes);
	running_stats->AddValues(STAT_PUSH_RESULTS_ADDED_BYTES, (int64_t) added_bytes);
}


// ************************************************************************************
// CMappingResultsDeliverer
// ************************************************************************************

// ************************************************************************************
CMappingResultsDeliverer::CMappingResultsDeliverer(CParams *params, CObjects *objects,  string _prefix_name, CJoinerMgr *_joiner_mgr)
{
	directory           = params->temp_prefix;
	prefix_name         = _prefix_name;
	q_res_ids			= objects->q_res_ids;
	joiner_mgr          = _joiner_mgr;

	keep_temporary_files = params->keep_temporary_files;

	mem_monitor         = objects->mem_monitor;
	serial_processing   = objects->serial_processing;

	running_stats       = objects->running_stats;

	verbosity_level     = params->verbosity_level;

	ptr_pool = objects->ptr_pool;
}

// ************************************************************************************
CMappingResultsDeliverer::~CMappingResultsDeliverer()
{
}

// ************************************************************************************
void CMappingResultsDeliverer::operator()()
{
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	uint32_t res_id;

	while(!q_res_ids->IsCompleted())
	{
		if(q_res_ids->Pop(res_id))
		{
			string file_name = get_name(res_id);

			shared_ptr<CMapperFile> res_group_file(new CMapperFile(EXT_MAPPING_GROUP, MARKER_MAPPING_GROUP));
			res_group_file->OpenRead(file_name);
			uint64_t file_size = res_group_file->GetSize();

			if(verbosity_level >= 2)
				cerr << "Reading res. group: " << res_id << " of size " << file_size << "\n";

			mem_monitor->ForceIncrease(file_size*2);

			uchar_t* data = ptr_pool->Allocate(file_size);
			serial_processing->Call([&]{res_group_file->Read(data, file_size); });
	
			res_group_file->Close();
			if (!keep_temporary_files)
				res_group_file->Remove();

			joiner_mgr->PutResultsBin(res_id, results_bin_t(data, file_size));
		}
	}

	if(verbosity_level >= 2)
		cerr << "MappingResultsDeliverer: Mark completed\n";
	joiner_mgr->MarkBinsCompleted();

	thr_watch.StopTimer();

	running_stats->AddTotals(STAT_TIME_THR_RES_READER, thr_watch.GetElapsedTime());
}

// ************************************************************************************
string CMappingResultsDeliverer::get_name(int n)
{
	return directory + prefix_name + Int2StringFilled(n, 5);
}


// ************************************************************************************
// CResultGroupsWriter
// ************************************************************************************

// ************************************************************************************
CResultGroupsWriter::CResultGroupsWriter(CParams *params, CObjects *objects, string _prefix_name)
{
	directory           = params->temp_prefix;
	prefix_name         = _prefix_name;
	mp_map_res          = objects->mp_map_res;
	q_map_res			= objects->q_map_res;
	no_groups		    = params->no_res_groups;
	min_no_free_parts   = params->min_no_free_group_parts;
	max_mem_single_file = params->max_group_part_size;

	running_stats       = objects->running_stats;

	writing_buffer	= nullptr;
	n_parts         = 0;

	SetMemcpyCacheLimit(8);

	// If alread registered nothing happens
	running_stats->Register(STAT_MAPPING_GROUP_WRITER_TOTAL, "Mapping groups writer: sizes (total)", running_stats_t::averages);
}

// ************************************************************************************
CResultGroupsWriter::~CResultGroupsWriter()
{
	if(writing_buffer)
		delete[] writing_buffer;
}

// ************************************************************************************
void CResultGroupsWriter::operator()()
{
	res_group_t group_part;
	
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	open_files();

	while(!q_map_res->IsCompleted())
	{
		if(q_map_res->Pop(group_part))
		{
			if(buffer_sizes[group_part.group_id] + group_part.size > max_mem_single_file)
				write_part(group_part.group_id);
			while(mp_map_res->GetAvailableParts() < (int64_t) min_no_free_parts && q_map_res->GetSize() < min_no_free_parts)
				write_part(largest_group_id);

			buffer_parts[group_part.group_id].push_back(group_part);
			buffer_sizes[group_part.group_id]  += group_part.size;

			if(buffer_sizes[group_part.group_id] > largest_group_size)
			{
				largest_group_size = buffer_sizes[group_part.group_id];
				largest_group_id   = group_part.group_id;
			}
			++n_parts;
		}
	}

	// Move all remaining parts to queue
	for(uint32_t i = 0; i < no_groups; ++i)
		write_part(i);

	close_files();

	thr_watch.StopTimer();

	running_stats->AddTotals(STAT_TIME_THR_RES_WRITER, thr_watch.GetElapsedTime());
}

// ************************************************************************************
// Return name of a file related to a kmer of given id. Also return disk no. for that kmer
string CResultGroupsWriter::get_name(int n)
{
	return directory + prefix_name + Int2StringFilled(n, 5);
}

// ************************************************************************************
bool CResultGroupsWriter::open_files()
{
	writing_buffer = new uchar_t[max_mem_single_file];

	f_groups.clear();

	for(uint32_t i = 0; i < no_groups; ++i)
	{
		f_groups.push_back(shared_ptr<CMapperFile>(new CMapperFile(EXT_MAPPING_GROUP, MARKER_MAPPING_GROUP)));
		f_groups.back()->OpenWrite(get_name(i));
	}

	buffer_sizes.resize(no_groups);
	fill(buffer_sizes.begin(), buffer_sizes.end(), 0);

	buffer_parts.resize(no_groups);
	for(auto &p : buffer_parts)
		p.clear();

	largest_group_size = 0;
	largest_group_id   = 0;
	n_parts            = 0;

	return true;
}

// ************************************************************************************
bool CResultGroupsWriter::close_files()
{
	f_groups.clear();
	delete[] writing_buffer;
	writing_buffer = nullptr;

	return true;
}

// ************************************************************************************
bool CResultGroupsWriter::write_part(uint32_t group_id)
{
	uint64_t offset = 0;

	// Copy all the parts to single memory place
	for(auto &p : buffer_parts[group_id])
	{
		A_memcpy(writing_buffer+offset, p.data, p.size);
		offset += p.size;
		mp_map_res->Free(p.data);
		--n_parts;
	}

	if(!offset)
		return true;

	// Write to disk
	f_groups[group_id]->Write(writing_buffer, offset);
	buffer_parts[group_id].clear();
	buffer_sizes[group_id] = 0;

	running_stats->AddValues(STAT_MAPPING_GROUP_WRITER_TOTAL, (int64_t) offset);

	// Find the largest bin
	largest_group_id   = 0;
	largest_group_size = buffer_sizes[0];
	for(uint32_t i = 1; i < no_groups; ++i)
		if(buffer_sizes[i] > largest_group_size)
		{
			largest_group_size = buffer_sizes[i];
			largest_group_id   = i;
		}

	return true;
}

// EOF
