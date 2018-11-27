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


#include "bins.h"
#include <algorithm>
#include <map>
#include <functional>
#include <queue>
#include <vector>
#include <string>
#include <cmath>
#include "../libs/asmlib.h"

using namespace std;

// ************************************************************************************
// CBinsWriter
// ************************************************************************************

// ************************************************************************************
CBinsWriter::CBinsWriter(CParams *params, CObjects *objects, string _name_prefix, uint32_t _stage_id, bool _final_substage)
{
	directory		  = params->temp_prefix;
	name_prefix		  = _name_prefix;
	mp_bins			  = objects->mp_bins_write;
	q_bins_write	  = objects->q_bins_write;
	running_stats     = objects->running_stats;
	stage_id		  = _stage_id;

	final_substage = _final_substage;

	prefix_file_names = params->prefix_file_names;

	no_bins           = (uint32_t) prefix_file_names.size();

	min_no_free_parts   = params->min_no_free_parts;
	max_mem_single_file = params->max_bin_part_size;

	keep_temporary_files = params->keep_temporary_files;

	writing_buffer = nullptr;

	n_total_parts = 0;
	n_write_bytes = 0;

	// If alread registered nothing happens
	running_stats->Register(STAT_BINS_WRITER_TOTAL, "Bin writer: sizes (total)", running_stats_t::averages);
	running_stats->Register(STAT_BINS_WRITER_BASE + stage_id, "Bin writer: sizes (" + StageDesc(stage_id) + ")", running_stats_t::averages);

	verbosity_level = params->verbosity_level;

	SetMemcpyCacheLimit(8);
}

// ************************************************************************************
CBinsWriter::~CBinsWriter()
{
	if(writing_buffer)
		delete[] writing_buffer;
}

// ************************************************************************************
void CBinsWriter::operator()()
{
	reads_bin_t bin_part;
	
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	open_files();

	while(!q_bins_write->IsCompleted())
	{
		if(q_bins_write->Pop(bin_part))
		{
			if(buffer_sizes[bin_part.bin_id] + bin_part.size > max_mem_single_file)
				write_bin(bin_part.bin_id);
			while(mp_bins->GetAvailableParts() < (int64_t) min_no_free_parts && q_bins_write->GetSize() < min_no_free_parts)
				write_bin(largest_bin_id);

			buffer_parts[bin_part.bin_id].push_back(bin_part);
			buffer_sizes[bin_part.bin_id]  += bin_part.size;
			read_counters[bin_part.bin_id] += bin_part.count;

			if(buffer_sizes[bin_part.bin_id] > largest_bin_size)
			{
				largest_bin_size = buffer_sizes[bin_part.bin_id];
				largest_bin_id   = bin_part.bin_id;
			}
			++n_parts;
			++n_total_parts;
		}
	}

	// Move all remaining parts to queue
	for(uint32_t i = 0; i < no_bins; ++i)
		write_bin(i);

	close_files();
	store_stats();

	thr_watch.StopTimer();

	running_stats->AddTotals(STAT_TIME_THR_BIN_WRITER_BASE + stage_id, thr_watch.GetElapsedTime());
}

// ************************************************************************************
// Return name of a file related to a kmer of given id. Also return disk no. for that kmer
string CBinsWriter::get_name(int n)
{
	return directory + name_prefix + prefix_file_names[n];
}

// ************************************************************************************
bool CBinsWriter::open_files()
{
	writing_buffer = new uchar_t[max_mem_single_file];

	f_bins.clear();

	if (!final_substage)
		for(uint32_t i = 0; i < no_bins; ++i)
		{
			f_bins.push_back(shared_ptr<CMapperFile>(new CMapperFile(EXT_READS_BIN, MARKER_READS_BIN)));
			f_bins.back()->OpenWrite(get_name(i));
		}

	buffer_sizes.resize(no_bins);
	fill(buffer_sizes.begin(), buffer_sizes.end(), 0);

	read_counters.resize(no_bins);
	fill(read_counters.begin(), read_counters.end(), 0);

	buffer_parts.resize(no_bins);
	for(auto &p : buffer_parts)
		p.clear();

	largest_bin_size = 0;
	largest_bin_id   = 0;
	n_parts			 = 0;

	return true;
}

// ************************************************************************************
void CBinsWriter::store_stats()
{
	if (!keep_temporary_files && final_substage)
		return;

	shared_ptr<CMapperFile> f_stats(new CMapperFile(EXT_BIN_STATS, MARKER_BIN_STATS));

	f_stats->OpenWrite(directory + name_prefix + "bin_stats");
	f_stats->Write(&no_bins, sizeof(uint32_t));
	for(auto &p: read_counters)
		f_stats->Write(&p, sizeof(uint64_t));

	f_stats->Close();
}

// ************************************************************************************
bool CBinsWriter::close_files()
{
	if (!final_substage)
		f_bins.clear();

	delete[] writing_buffer;
	writing_buffer = nullptr;

	return true;
}

// ************************************************************************************
bool CBinsWriter::write_bin(uint32_t bin_id)
{
	uint64_t offset = 0;

	// Copy all the parts to single memory place
	for(auto &p : buffer_parts[bin_id])
	{
		A_memcpy(writing_buffer+offset, p.data, p.size);
		offset += p.size;
		mp_bins->Free(p.data);
		--n_parts;
	}

	if(!offset)
		return true;

	// Write to disk
	if(!final_substage)
		f_bins[bin_id]->Write(writing_buffer, offset);
	buffer_parts[bin_id].clear();
	buffer_sizes[bin_id] = 0;

	running_stats->AddValues(STAT_BINS_WRITER_TOTAL, (int64_t) offset);
	running_stats->AddValues(STAT_BINS_WRITER_BASE + stage_id, (int64_t) offset);

	n_write_bytes += offset;
	if(verbosity_level >= 3)
		cerr << "Write bin: " << offset << " B\n";

	// Find the largest bin
	largest_bin_id   = 0;
	largest_bin_size = buffer_sizes[0];
	for(uint32_t i = 1; i < no_bins; ++i)
		if(buffer_sizes[i] > largest_bin_size)
		{
			largest_bin_size = buffer_sizes[i];
			largest_bin_id   = i;
		}

	return true;
}


// ************************************************************************************
// CBinsReader
// ************************************************************************************

// ************************************************************************************
CBinsReader::CBinsReader(CParams *params, CObjects *objects, string _name_prefix, uint32_t _stage_id)
{
	directory             = params->temp_prefix;
	name_prefix           = _name_prefix;
	prefix_file_names     = params->prefix_file_names;

	running_stats         = objects->running_stats;

	mem_monitor		      = objects->mem_monitor;
	q_bins_read           = objects->q_bins_read;

	stage_id              = _stage_id;

	max_reads_compression = params->max_reads_compression;
	read_len		      = (uint32_t) params->read_len;

	verbosity_level       = params->verbosity_level;

	keep_temporary_files = params->keep_temporary_files;
}

// ************************************************************************************
CBinsReader::~CBinsReader()
{
}

// ************************************************************************************
void CBinsReader::operator()()
{
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	if(!load_stats())
	{
		cerr << "Error: Cannot load bin stats\n";
		exit(1);
	}

	shared_ptr<CMapperFile> f_bin(new CMapperFile(EXT_READS_BIN, MARKER_READS_BIN));

	// Change order of bins according to decreasing no. of reads
	vector<pair<uint32_t, size_t>> bin_reordering;
	for(uint32_t i = 0; i < no_bins; ++i)
		bin_reordering.push_back(make_pair(i, read_counters[i]));
	
	// Sort bins from the largest one to the smallest, but put mono-letter bins first
	sort(bin_reordering.begin(), bin_reordering.end(), [&](pair<uint32_t, size_t> x, pair<uint32_t, size_t> y) {
		size_t x_cnt = x.second;
		size_t y_cnt = y.second;

		if (is_mono_symbol(x.first))
			x_cnt <<= 32;
		if (is_mono_symbol(y.first))
			y_cnt <<= 32;

		return x_cnt > y_cnt;
	});

	for(uint32_t i = 0; i < no_bins; ++i)
	{
		uint32_t bin_id = bin_reordering[i].first;
		if(verbosity_level >= 3)		
			cerr << "Read bin: " << get_name(bin_id) << "  (" << i << ")\n";
		if(!f_bin->OpenRead(get_name(bin_id)))
			continue;

		uint64_t raw_size = f_bin->GetSize();
		uint64_t size;
		uint64_t count = read_counters[bin_id];
		uint64_t req_memory;

		// Wait for possibility of allocation of necessary amount of memory
		if(max_reads_compression)
			// More memory will be necessary since reads decompression will be made
			size = (raw_size * 12 + 6) / 7 + read_len + 32;
		else
			size = raw_size;
		req_memory = size + count * sizeof(uint32_t);

		mem_monitor->ForceIncrease(req_memory);
		uchar_t *data = new uchar_t[size];
		uchar_t *raw_data = data + size - raw_size;

		f_bin->Read(raw_data, raw_size);
		q_bins_read->Push(reads_bin_t(i, bin_id, data, raw_size, size, read_counters[bin_id]));

		f_bin->Close();

		if(!keep_temporary_files)
			f_bin->Remove();
	}

	q_bins_read->MarkCompleted();

	thr_watch.StopTimer();

	running_stats->AddTotals(STAT_TIME_THR_BIN_READER_BASE + stage_id, thr_watch.GetElapsedTime());
}

// ************************************************************************************
string CBinsReader::get_name(int n)
{
	return directory + name_prefix + prefix_file_names[n];
}

// ************************************************************************************
bool CBinsReader::is_mono_symbol(int n)
{
	string &s = prefix_file_names[n];

	int cnt = 0;

	for (char c : {'A', 'C', 'G', 'T'})
		cnt += find(s.begin(), s.end(), c) != s.end();

	return cnt == 1;
}

// ************************************************************************************
bool CBinsReader::load_stats()
{
	shared_ptr<CMapperFile> f_stats(new CMapperFile(EXT_BIN_STATS, MARKER_BIN_STATS));

	if(!f_stats->OpenRead(directory + name_prefix + "bin_stats"))
		return false;

	f_stats->Read(&no_bins, sizeof(uint32_t));
	read_counters.resize(no_bins);	
	for(auto &p: read_counters)
		f_stats->Read(&p, sizeof(uint64_t));

	f_stats->Close();

	if (!keep_temporary_files)
		f_stats->Remove();

	return true;
}


// ************************************************************************************
// CBinPrefixes
// ************************************************************************************

// ************************************************************************************
//
CBinPrefixes::CBinPrefixes(uint32_t _max_prefix_len, uint32_t _verbosity_level)
{
	max_prefix_len = _max_prefix_len;
	prefix_table_size = (1 << (2 * max_prefix_len));
	is_valid = false;

	verbosity_level = _verbosity_level;
}
	
// ************************************************************************************
CBinPrefixes::~CBinPrefixes()
{
}

// ************************************************************************************
// Set counts of prefixes using LUT from SA
void CBinPrefixes::SetLUT(uint32_t *lut_data, uint32_t lut_prefix_len)
{
	// Count occurrences of prefixes of length max_prefix_len
	small_lut.resize(prefix_table_size+1);
	prefix_map.resize(prefix_table_size);

	uint32_t fill_shift = 2*lut_prefix_len - 2*max_prefix_len;
	for(int i = 0; i <= (1 << (2*max_prefix_len)); ++i)
		small_lut[i] = lut_data[i << fill_shift];
}

// ************************************************************************************
// Reduces no. of bins to the given value
bool CBinPrefixes::ReduceMap(uint32_t _no_bins)
{
	uint32_t current_no_bins = 1 << (2 * max_prefix_len);

	if (small_lut.empty())
		return false;

	decodes['A'] = 0;		codes[0] = 'A';
	decodes['C'] = 1;		codes[1] = 'C';
	decodes['G'] = 2;		codes[2] = 'G';
	decodes['T'] = 3;		codes[3] = 'T';

	no_bins = _no_bins;

	typedef pair<string, uint64_t> pfx_desc_t;
	auto cmp = [](pfx_desc_t x, pfx_desc_t y) {return x.second > y.second; };
	
	// pq stores prefixes with sizes - smallest first
	priority_queue<pfx_desc_t, vector<pfx_desc_t>, decltype(cmp)> pq(cmp);		

	// bins for which not all strings are 
	map<string, uint64_t> pfx_active;

	// Construct map of prefix strings related with counters
	uint64_t sum = 0;
	for (uint32_t i = 0; i < current_no_bins; ++i)
	{
		// decode prefix string and its rev. compl.
		string prefix = "";
		uint32_t x = i;

		for (uint32_t j = 0; j < max_prefix_len; ++j)
		{
			prefix.push_back(codes[x & 3]);
			x >>= 2;
		}
		reverse(prefix.begin(), prefix.end());

		uint64_t cnt = count_suffixes_from(prefix);

		cnt = (uint64_t) pow((double) cnt, 1.5);
		pfx_active[prefix] = cnt;

		sum += cnt;
		if (i % 4 == 3)
		{
			string str = prefix;
			str.pop_back();

			pq.push(make_pair(str, sum));
			sum = 0;
		}
	}
	
	// Reduce no. of bins
	while (current_no_bins >= no_bins)
	{
		pfx_desc_t x = pq.top();
		pq.pop();
		string str = x.first;

		for (string c : {"A", "C", "G", "T"})
			pfx_active.erase(str + c);

		pfx_active[str] = x.second;
		current_no_bins -= 3;

		if (str.size() == 3)		// do not allow to create shorter than 3-symbol long prefixes
			continue;

		str.pop_back();
		int cnt_active = 0;
		for (string c : {"A", "C", "G", "T"})
			cnt_active += (int) pfx_active.count(str + c);

		if (cnt_active == 4)
		{
			uint64_t sum_str = 0;
			for (string c : {"A", "C", "G", "T"})
				sum_str += pfx_active[str + c];

			pq.push(make_pair(str, sum_str));
		}

		if (verbosity_level >= 3)
			cerr << "Bin no. reducer: " << current_no_bins << "\n";
	}

	// Construct map of prefixes onto bin ids
	prefix_file_names.clear();

	uint32_t i = 0;
	for (auto p = pfx_active.begin(); p != pfx_active.end(); ++p, ++i)
	{
		prefix_file_names.push_back(p->first);
		uint32_t p_len = (uint32_t)p->first.length();

		uint32_t pref_decoded = 0;
		for (uint32_t j = 0; j < p_len; ++j)
			pref_decoded = (pref_decoded << 2) + decodes[(uchar_t)p->first[j]];

		uint32_t fill_shift = 2 * (max_prefix_len - p_len);
		for (uint32_t j = pref_decoded << fill_shift; j < (pref_decoded + 1) << fill_shift; ++j)
			prefix_map[j] = i;
	}
	prefix_file_names.push_back("X");

	is_valid = true;

	return true;
}

// ************************************************************************************
// Compute the no. of suffixes starting from (includes both direct and rev. compl.)
uint64_t CBinPrefixes::count_suffixes_from(string str)
{
	uint64_t res;
	uint32_t shift_to_fill = 2*(max_prefix_len - (uint32_t) str.length());

	// Decode suffix
	uint32_t dir_bin = 0;
	uint32_t rc_bin  = 0;
	for(uint32_t i = 0; i < str.length(); ++i)
	{
		dir_bin = (dir_bin << 2) + decodes[(uchar_t) str[i]];
		rc_bin  = rc_bin + ((3-decodes[(uchar_t) str[i]]) << (2*i));
	}

	res  = small_lut[(dir_bin+1) << shift_to_fill] - small_lut[dir_bin << shift_to_fill];		// direct suffixes
	res += small_lut[(rc_bin +1) << shift_to_fill] - small_lut[rc_bin  << shift_to_fill];		// rev. compl. suffixes

	return res;
}

// ************************************************************************************
// Return the map of prefixes and the file names of bins
bool CBinPrefixes::GetMaps(vector<uint32_t> &_prefix_map, vector<string> &_prefix_file_names)
{
	if(!is_valid)
		return false;

	_prefix_map.assign(prefix_map.begin(), prefix_map.end());
	_prefix_file_names.assign(prefix_file_names.begin(), prefix_file_names.end());

	return true;	
}

// EOF
