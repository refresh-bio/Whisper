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
// candidate_mapping_t
// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::lev>
	(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::lev;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	
	t.penalty = t.calc_penalty_lev();
	
	return t;
}

// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::clipping_mismatches>
	(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, uint32_t _clipping_len)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::clipping_mismatches;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	t.clipping_len = _clipping_len;

	t.penalty = t.calc_penalty_mismatches_clipping();

	return t;
}

// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::mismatches_clipping>
(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, uint32_t _clipping_len)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::mismatches_clipping;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	t.clipping_len = _clipping_len;

	t.penalty = t.calc_penalty_mismatches_clipping();

	return t;
}

// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::indel_clipping>
(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, uint32_t _matching_len1,
	int32_t _indel_len1, uint32_t _matching_len2)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::indel_clipping;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	t.matching_len1 = _matching_len1;
	t.indel_len1 = _indel_len1;
	t.matching_len2 = _matching_len2;

	t.penalty = t.calc_penalty_indel_clipping();

	return t;
}

// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::clipping_indel>
(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, uint32_t _matching_len1,
	int32_t _indel_len1, uint32_t _matching_len2)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::clipping_indel;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	t.matching_len1 = _matching_len1;
	t.indel_len1 = _indel_len1;
	t.matching_len2 = _matching_len2;

	t.penalty = t.calc_penalty_indel_clipping();

	return t;
}

// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::indel1>
(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, uint32_t _matching_len1,
	int32_t _indel_len1)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::indel1;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	t.matching_len1 = _matching_len1;
	t.indel_len1 = _indel_len1;

	t.penalty = t.calc_penalty_indel1();
	
	return t;
}

// ************************************************************************************
template<>
candidate_mapping_t candidate_mapping_t::construct<mapping_type_t::indel2>
(ref_pos_t _pos, genome_t _direction, uint32_t _no_mismatches, uint32_t _matching_len1,
	int32_t _indel_len1, uint32_t _matching_len2, int32_t _indel_len2)
{
	candidate_mapping_t t;

	t.type = mapping_type_t::indel2;
	t.pos = _pos;
	t.direction = _direction;
	t.no_mismatches = _no_mismatches;
	t.matching_len1 = _matching_len1;
	t.indel_len1 = _indel_len1;
	t.matching_len2 = _matching_len2;
	t.indel_len2 = _indel_len2;

	t.penalty = t.calc_penalty_indel2();

	return t;
}

// ************************************************************************************
size_t candidate_mapping_t::serialize(uchar_t* ptr)
{
	auto ptr0 = ptr;

	StoreUInt(ptr, pos, 4);
	ptr += 4;

	*ptr = (direction == genome_t::direct) ? 0x80u : 0;
	*ptr += (type == mapping_type_t::lev) ? 0x40u : 0;
	*ptr++ += (uchar_t)(no_mismatches & 0x3f);
	   
	if (type == mapping_type_t::lev)
		return ptr - ptr0;

	*ptr++ = (uchar_t)type;

	if (type == mapping_type_t::indel1)
	{
		*ptr++ = (uchar_t) ((matching_len1 >> 8) & 0xf);
		*ptr++ = (uchar_t) (matching_len1 & 0xff);
		*ptr++ = indel_size_encode(indel_len1);
 	}
	else if (type == mapping_type_t::indel2)
	{
		*ptr = (uchar_t)((matching_len1 >> 8) & 0xf);
		*ptr++ += (uchar_t)((matching_len2 >> 4) & 0xf0u);

		*ptr++ = (uchar_t)(matching_len1 & 0xffu);
		*ptr++ = (uchar_t)(matching_len2 & 0xffu);
		*ptr++ = indel_size_encode(indel_len1);
		*ptr++ = indel_size_encode(indel_len2);
	}
	else if (type == mapping_type_t::indel_clipping || type == mapping_type_t::clipping_indel)
	{
		*ptr = (uchar_t)((matching_len1 >> 8) & 0xf);
		*ptr++ += (uchar_t)((clipping_len >> 4) & 0xf0u);

		*ptr++ = (uchar_t)(matching_len1 & 0xffu);
		*ptr++ = (uchar_t)(clipping_len & 0xffu);
		*ptr++ = indel_size_encode(indel_len1);
	}
	else if (type == mapping_type_t::mismatches_clipping || type == mapping_type_t::clipping_mismatches)
	{
		*ptr++ = (uchar_t)((matching_len1 >> 8) & 0xf);
		*ptr++ = (uchar_t)(matching_len1 & 0xffu);
	}
	else if (type == mapping_type_t::clipping_clipping)
	{
		*ptr = (uchar_t)((matching_len1 >> 8) & 0xf);
		*ptr++ += (uchar_t)((clipping_len >> 4) & 0xf0u);

		*ptr++ = (uchar_t)(matching_len1 & 0xffu);
		*ptr++ = (uchar_t)(clipping_len & 0xffu);
	}

	return ptr - ptr0;
}

// ************************************************************************************
size_t candidate_mapping_t::deserialize(uchar_t* ptr)
{
	auto ptr0 = ptr;

	uint64_t tmp;

	LoadUInt(ptr, tmp, 4);
	pos = (ref_pos_t)tmp;
	ptr += 4;

	direction = (*ptr & 0x80u) ? genome_t::direct : genome_t::rev_comp;
	no_mismatches = *ptr & 0x3f;
	if (*ptr++ & 0x40u)
	{
		type = mapping_type_t::lev;

		return ptr - ptr0;
	}

	type = (mapping_type_t)*ptr++;

	if (type == mapping_type_t::indel1)
	{
		matching_len1 = ((uint32_t) *ptr++) << 8;
		matching_len1 += (uint32_t) *ptr++;
		indel_len1 = indel_size_decode(*ptr++);
	}
	else if (type == mapping_type_t::indel2)
	{
		matching_len1 = ((uint32_t)(*ptr & 0xf)) << 8;
		matching_len2 = ((uint32_t)(*ptr++ & 0xf0u)) << 4;
		matching_len1 += (uint32_t)*ptr++;
		matching_len2 += (uint32_t)*ptr++;
		indel_len1 = indel_size_decode(*ptr++);
		indel_len2 = indel_size_decode(*ptr++);
	}
	else if (type == mapping_type_t::indel_clipping || type == mapping_type_t::clipping_indel)
	{
		matching_len1 = ((uint32_t)(*ptr & 0xf)) << 8;
		clipping_len = ((uint32_t)(*ptr++ & 0xf0)) << 4;
		matching_len1 += (uint32_t)*ptr++;
		clipping_len += (uint32_t)*ptr++;
		indel_len1 = indel_size_decode(*ptr++);
	}
	else if (type == mapping_type_t::mismatches_clipping || type == mapping_type_t::clipping_mismatches)
	{
		matching_len1 = ((uint32_t)(*ptr++ & 0xf)) << 8;
		matching_len1 += (uint32_t)*ptr++;
	}
	else if (type == mapping_type_t::clipping_clipping)
	{
		matching_len1 = ((uint32_t)(*ptr & 0xf)) << 8;
		clipping_len = ((uint32_t)(*ptr++ & 0xf0)) << 4;
		matching_len1 += (uint32_t)*ptr++;
		clipping_len += (uint32_t)*ptr++;
	}

	calc_penalty();

	return ptr - ptr0;
}

// ************************************************************************************
size_t candidate_mapping_t::check_rec_size(uchar_t* ptr)
{
	if (*(ptr + 4) & 0x40)		// lev
		return 5;

	mapping_type_t type = (mapping_type_t)*(ptr + 5);

	if (type == mapping_type_t::indel1)
		return 9;
	else if (type == mapping_type_t::indel2)
		return 11;
	else if (type == mapping_type_t::indel_clipping || type == mapping_type_t::clipping_indel)
		return 10;
	else if (type == mapping_type_t::mismatches_clipping || type == mapping_type_t::clipping_mismatches)
		return 8;
	else if (type == mapping_type_t::clipping_clipping)
		return 9;

	return 0;	// Cannot be here
}

// ************************************************************************************
uchar_t* candidate_mapping_t::skip(uchar_t* ptr)
{
	return ptr + check_rec_size(ptr);
}

// ************************************************************************************
void candidate_mapping_t::construct(candidate_mapping_t& mapping)
{
	type = mapping.type;
	pos = mapping.pos;
	direction = mapping.direction;
	no_mismatches = mapping.no_mismatches;
	matching_len1 = mapping.matching_len1;
	matching_len2 = mapping.matching_len2;
	indel_len1 = mapping.indel_len1;
	indel_len2 = mapping.indel_len2;
	clipping_len = mapping.clipping_len;

	calc_penalty();
}

// ************************************************************************************
void candidate_mapping_t::calc_penalty()
{
	if (type == mapping_type_t::lev)
		penalty = calc_penalty_lev();
	else if (type == mapping_type_t::clipping_indel || type == mapping_type_t::indel_clipping)
		penalty = calc_penalty_indel_clipping();
	else if (type == mapping_type_t::clipping_mismatches || type == mapping_type_t::mismatches_clipping)
		penalty = calc_penalty_mismatches_clipping();
	else if (type == mapping_type_t::indel1)
		penalty = calc_penalty_indel1();
	else if (type == mapping_type_t::indel2)
		penalty = calc_penalty_indel2();
	else if (type == mapping_type_t::clipping_clipping)
		penalty = calc_penalty_clipping_clipping();
	else
		penalty = 0;		// !! Should never be here
}

// ************************************************************************************
candidate_mapping_t& candidate_mapping_t::operator=(const candidate_mapping_t& mapping)
{
	if (&mapping == this)
		return *this;

	type = mapping.type;
	pos = mapping.pos;
	direction = mapping.direction;
	no_mismatches = mapping.no_mismatches;
	matching_len1 = mapping.matching_len1;
	matching_len2 = mapping.matching_len2;
	indel_len1 = mapping.indel_len1;
	indel_len2 = mapping.indel_len2;
	clipping_len = mapping.clipping_len;

	penalty = mapping.penalty;

	return *this;
}

// ************************************************************************************
bool candidate_mapping_t::operator<(const candidate_mapping_t& mapping)
{
	if (penalty != mapping.penalty)
		return penalty < mapping.penalty;

	return pos < mapping.pos;
}

// ************************************************************************************
candidate_mapping_t candidate_mapping_t::better(candidate_mapping_t& x, candidate_mapping_t& y)
{
	if (x.type == mapping_type_t::none)
		return y;
	if (y.type == mapping_type_t::none)
		return x;

	if (x < y)
		return x;
	else
		return y;
}

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
//void CMappingsHeapGatherer::Push(uint32_t pos, genome_t direction, uint32_t no_errors, int32_t indel_marker)
void CMappingsHeapGatherer::Push(candidate_mapping_t &mapping)
{
	if (size < max_size)
	{
		v_mappings[size++].construct(mapping);
		push_heap(v_mappings.begin(), v_mappings.begin() + size);
	}
	else
	{
		if (mapping < v_mappings.front())
		{
			pop_heap(v_mappings.begin(), v_mappings.end());
			v_mappings.back().construct(mapping);
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
candidate_mapping_t CMappingsHeapGatherer::Pop()
{
	candidate_mapping_t r;
		
	if (size)
	{
		r = v_mappings.front();
		pop_heap(v_mappings.begin(), v_mappings.begin() + size);
		--size;
	}

	return r;
}

// ************************************************************************************
candidate_mapping_t CMappingsHeapGatherer::PopUnsorted()
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
	return (uint32_t)(x >> 40) & 0xff;
}

// ************************************************************************************
uint32_t CMappingsHeapGatherer::DecodeIndelMarker(uint64_t x) const
{
	return (uint32_t)(x >> 48);
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

	buffer_reserve = id_bytes + mapping_counter_size + candidate_mapping_t::max_size;

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
void CMappingResultsCollector::Push(read_id_t id, candidate_mapping_t mapping)
	//ref_pos_t pos, genome_t direction, uint32_t no_differences, int32_t indel_marker)
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
		cur_read_counter_ptr = buffers[group_id] + buffer_sizes[group_id];

		buffer_sizes[group_id] += mapping_counter_size;
		added_bytes += mapping_counter_size;
		prev_read_id    = id;
		cur_read_count  = 0;
		true_read_count = 0;

		push_unique_counter++;
	}

//	IncrementUInt(buffers[group_id]+(buffer_sizes[group_id] - mapping_counter_size - cur_read_count*rec_size), mapping_counter_size);
	IncrementUInt(cur_read_counter_ptr, mapping_counter_size);

	auto serialized_size = mapping.serialize(buffers[group_id] + buffer_sizes[group_id]);
	buffer_sizes[group_id] += serialized_size;
	added_bytes += serialized_size;
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

#ifdef COLLECT_STATS
	running_stats->AddValues(STAT_PUSH_RESULTS    , (int64_t) push_counter);
	running_stats->AddValues(STAT_PUSH_RESULTS_UNQ, (int64_t) push_unique_counter);
	running_stats->AddValues(STAT_PUSH_RESULTS_SEND_BYTES, (int64_t) sent_bytes);
	running_stats->AddValues(STAT_PUSH_RESULTS_ADDED_BYTES, (int64_t) added_bytes);
#endif
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

#ifdef COLLECT_STATS
	running_stats->AddTotals(STAT_TIME_THR_RES_READER, thr_watch.GetElapsedTime());
#endif
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
#ifdef COLLECT_STATS
	running_stats->Register(STAT_MAPPING_GROUP_WRITER_TOTAL, "Mapping groups writer: sizes (total)", running_stats_t::averages);
#endif
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

#ifdef COLLECT_STATS
	running_stats->AddTotals(STAT_TIME_THR_RES_WRITER, thr_watch.GetElapsedTime());
#endif
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

#ifdef COLLECT_STATS
	running_stats->AddValues(STAT_MAPPING_GROUP_WRITER_TOTAL, (int64_t) offset);
#endif

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
