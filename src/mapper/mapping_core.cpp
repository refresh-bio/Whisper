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


#include "mapping_core.h"
#include "vector_utils.h"

#include "LevMyers128.h"
#include "LevMyers256.h"

#include <algorithm>
#include <utility>
#include <assert.h>

//**********************************************************************************************************
//
//**********************************************************************************************************
CMappingCore::CMappingCore(CParams *params, CObjects *objects, uint32_t _stage_id, mapping_mode_t _mapping_mode)
{
	instruction_set_t instruction_set;
	int x = instrset_detect();
	if (x >= 0 && x <= 8)
		instruction_set = (instruction_set_t)x;
	else if (x < 0)
		instruction_set = instruction_set_t::none;
	else
		instruction_set = instruction_set_t::avx2;
	
	prefix_map = params->prefix_map;
	q_bins_read = objects->q_bins_read;
	q_bins_write = objects->q_bins_write;
	q_map_res = objects->q_map_res;

	running_stats = objects->running_stats;
	stage_id = _stage_id;
	mapping_mode = _mapping_mode;

	mem_monitor = objects->mem_monitor;
	mp_bins_write = objects->mp_bins_write;
	mp_map_res = objects->mp_map_res;

	verbosity_level = params->verbosity_level;

	stage_major = (_stage_id >> 5) & 0x1f;
	stage_minor = _stage_id & 0x1f;
	sensitive_mode = (bool)(_stage_id >> 10);
	sensitivity_factor = params->sensitivity_factor;
	max_no_mappings = params->max_no_mappings;

	no_bins = prefix_map.back() + 2;		// the last element of prefix map contains the last bin_id
	read_id_len = BITS2BYTES(bin_total_bits);
	start_pos = params->stage_segments[stage_major][stage_minor].first;
	end_pos = params->stage_segments[stage_major][stage_minor].second;

	if (stage_minor < stage_major)
	{
		next_start_pos = params->stage_segments[stage_major][stage_minor + 1].first;
		next_end_pos = params->stage_segments[stage_major][stage_minor + 1].second;
	}
	else
	{
		next_start_pos = params->stage_segments[stage_major + 1][0].first;
		next_end_pos = params->stage_segments[stage_major + 1][0].second;
	}
	params_aux_ptr = params;

	for (uint32_t l = 0; l < stage_minor; l++)
		segment_segments.push_back(make_pair(params->stage_segments[stage_major][l].first, params->stage_segments[stage_major][l].second));

	no_groups = params->no_res_groups;

	prefixes = params->prefix_file_names;
	sa_prefix_overhead = params->sa_prefix_overhead;

	bin_prefix = IntLog4(prefix_map.size());

	reference = objects->reference;
	sa_dir = objects->sa_dir;
	sa_rc = objects->sa_rc;
	sa_in_ram = params->sa_in_ram;

	max_read_len = params->max_read_len;
	min_read_len = params->min_read_len;

	old_sequence = new uchar_t[max_read_len / 2 + max_read_len % 2];
	old_sequence_sft = new uchar_t[max_read_len / 2 + 1];

	CUR_OFFSET = start_pos;
	predef_match_length = end_pos - start_pos;

	no_of_substrings = stage_major + 5;

	// Assure that the substrings are not longer (in total) than reads 
	while (no_of_substrings * max_substr_len > min_read_len)
		--no_of_substrings;

	if (no_of_substrings < stage_major + 1)
		no_of_substrings = stage_major + 1;

#ifdef USE_128b_SSE_CF_AND_DP			
	uint32_t max_no_of_substrings_for_DP = 0;

	for (uint32_t i = 0; i < stage_major; i++)	//stage_major corresponds to the number of mismatches
		max_no_of_substrings_for_DP += (2 * i + 1);

	max_no_of_substrings_for_DP += (no_of_substrings - stage_major)*(2 * stage_major + 1);

	dir_aux_arr = new CAux_arrays[max_no_of_substrings_for_DP];
	rc_aux_arr = new CAux_arrays[max_no_of_substrings_for_DP];
#else
	dir_aux_arr = new CAux_arrays[no_of_substrings];
	rc_aux_arr = new CAux_arrays[no_of_substrings];
#endif
	dir_aux_params = new CAux_params;
	rc_aux_params = new CAux_params;

	adjust_substr_pos_arr();

	cur_substr_values_from_pattern = new uint32_t[no_of_substrings];

	reads_collector = new CReadsCollector(params, objects);
	reads_deliverer = new CReadsDeliverer(params, read_id_len, start_pos, end_pos, running_stats, objects->mp_reads, stage_id);
	results_collector = new CMappingResultsCollector(params, objects, stage_id);

	progress = objects->progress;

#ifdef USE_128b_SSE_CF_AND_DP
	prev_dp_ptr = new uint32_t[max_read_len + 1 + 2 * stage_major + 2];
	curr_dp_ptr = new uint32_t[max_read_len + 1 + 2 * stage_major + 2];
	genome_prefetch = new uchar_t[max_read_len + 1 + 2 * stage_major + 2];

	edit_dist = new CEditDist((uint32_t)max_read_len, (uint32_t)max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2,
		(uint32_t)(stage_major*params->sensitivity_factor));
	edit_dist->SetReference(reference->GetData(), sa_dir->GetRefSize(), start_pos);

	levMyers64 = new LevMyers64((uint32_t) max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2), 0);

	switch (instruction_set) {
	case instruction_set_t::none:
	case instruction_set_t::sse:
		throw new std::runtime_error("SSE2 extensions required!");
	case instruction_set_t::sse2:
	case instruction_set_t::sse3:
	case instruction_set_t::sse3s:
	case instruction_set_t::sse41:
	case instruction_set_t::sse42:
		levMyers128 = new LevMyers128<instruction_set_t::sse2>((uint32_t) max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2),	0);
		levMyers256 = new LevMyers64((uint32_t)max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2), 0);
		break;
	case instruction_set_t::avx:
		levMyers128 = new LevMyers128<instruction_set_t::avx>((uint32_t) max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2), 0);
		levMyers256 = new LevMyers64((uint32_t) max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2), 0);
		break;
	case instruction_set_t::avx2:
		levMyers128 = new LevMyers128<instruction_set_t::avx2>((uint32_t) max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2),	0);
		levMyers256 = new LevMyers256((uint32_t) max_read_len, (uint32_t) (max_read_len + 2 * (uint32_t)(stage_major*params->sensitivity_factor) + 2), 0);
		break;
	}

	levMyers64->setReference(reference->GetData(), sa_dir->GetRefSize(), start_pos);
	levMyers128->setReference(reference->GetData(), sa_dir->GetRefSize(), start_pos);
	levMyers256->setReference(reference->GetData(), sa_dir->GetRefSize(), start_pos);

	CF_array_dir_gen = new uint32_t[1];
	CF_array_rc_gen = new uint32_t[1];

	for (int32_t i = 0; i < 256; ++i)
		CF_SSE_lut[i] = Vec16c(0);

	CF_SSE_lut[0x00].load(list[0]);			// AA
	CF_SSE_lut[0x01].load(list[1]);			// AC
	CF_SSE_lut[0x02].load(list[2]);			// AG
	CF_SSE_lut[0x03].load(list[3]);			// AT

	CF_SSE_lut[0x10].load(list[4]);			// CA
	CF_SSE_lut[0x11].load(list[5]);			// CC
	CF_SSE_lut[0x12].load(list[6]);			// CG
	CF_SSE_lut[0x13].load(list[7]);			// CT

	CF_SSE_lut[0x20].load(list[8]);			// GA
	CF_SSE_lut[0x21].load(list[9]);			// GC
	CF_SSE_lut[0x22].load(list[10]);		// GG
	CF_SSE_lut[0x23].load(list[11]);		// GT

	CF_SSE_lut[0x30].load(list[12]);		// TA
	CF_SSE_lut[0x31].load(list[13]);		// TC
	CF_SSE_lut[0x32].load(list[14]);		// TG
	CF_SSE_lut[0x33].load(list[15]);		// TT
#endif
}

//**********************************************************************************************************
//
//**********************************************************************************************************
CMappingCore::~CMappingCore()
{
	delete reads_collector;
	delete reads_deliverer;
	delete results_collector;

	delete[] dir_aux_arr;
	delete[] rc_aux_arr;
	delete dir_aux_params;
	delete rc_aux_params;
	delete[] old_sequence;
	delete[] old_sequence_sft;
	delete[] cur_substr_values_from_pattern;

#ifdef USE_128b_SSE_CF_AND_DP
	delete[] CF_array_dir_gen;
	delete[] CF_array_rc_gen;

	delete[] prev_dp_ptr;
	delete[] curr_dp_ptr;
	delete[] genome_prefetch;

	delete edit_dist;

	delete levMyers256;
	delete levMyers128;
	delete levMyers64;
#endif
}

//**********************************************************************************************************
//
//**********************************************************************************************************
void CMappingCore::operator()()
{
	n_reads = 0;

	CThreadWatch thr_watch, bin_watch;

	thr_watch.StartTimer();

	// Mapping reads
	while (!q_bins_read->IsCompleted())
	{
		reads_bin_t bin;
		if (q_bins_read->Pop(bin))
		{
			bin_watch.StartTimer();
			
			if (verbosity_level >= 2)
			{
				cerr << "Mapping bin: " << prefixes[bin.bin_id] << "   (" << bin.bin_num << ": " << bin.count << " reads)\n";
			}

			prefix_len = (uint32_t)prefixes[bin.bin_id].size();

			if (bin.bin_id < no_bins - 1)
				adjust_group_sizes(bin);

			reads_deliverer->SetBin(bin);
			reads_deliverer->Start();
			if (bin.bin_id < no_bins - 1)
				if ((stage_major == 1) && !stage_minor)
					process_reads_init_stage(bin);
				else
					process_reads(bin);
			else
				copy_reads(bin);		// do not map reads with Ns in prefix

			delete[] bin.data;
			mem_monitor->Decrease(bin.size + bin.count * sizeof(uint32_t));

			progress->Step(1);

			if (verbosity_level >= 2)
			{
				cerr << "Completed bin: " << prefixes[bin.bin_id] << "\n";
			}
			bin_watch.StopTimer();

			uint32_t stat_id = STAT_BIN_TIMES + stage_id * 1000 + bin.bin_id;
			running_stats->Register(stat_id, "Bin times " + StageDesc(stage_id) + " : " + prefixes[bin.bin_id] + " : ", running_stats_t::totals);
			running_stats->AddTotals(stat_id, bin_watch.GetElapsedTime());
		}
	}

	complete();

	q_bins_write->MarkCompleted();
	q_map_res->MarkCompleted();

	thr_watch.StopTimer();

	running_stats->AddValues(STAT_TIME_THR_MAPPING_CORE_BASE + stage_id, thr_watch.GetElapsedTime());
}

//**********************************************************************************************************
//
//**********************************************************************************************************
inline uint32_t CMappingCore::get_next_bin_id(uchar_t* data)
{
	uint32_t bin_id = 0;
	uint32_t i;

	// Look for Ns
	i = next_start_pos;
	if (i & 1)		// even starting position in a read
	{
		if (LO_NIBBLE(data[i / 2]) & 0x44)
			return no_bins - 1;
		i++;
	}
	for (; i + 1 < next_end_pos; i += 2)
		if (data[i / 2] & 0x44)
			return no_bins - 1;

	if (next_end_pos & 1)
		if (HI_NIBBLE(data[i / 2] & 0x44))
			return no_bins - 1;

	// Compute bin id
	for (i = 0; i < bin_prefix; ++i)
		bin_id = (bin_id << 2) + GET_SYMBOL(data, next_start_pos + i);

	bin_id = prefix_map[bin_id];		// maps long prefix of a read onto bin id

	return bin_id;
}

//**********************************************************************************************************
// Return binary encoded (2bits per symbol) size-symbol long prefix of read
//********************************************************************************************************** 
inline uint32_t CMappingCore::get_read_prefix(uchar_t* data, uint32_t size)
{
	uint32_t r = 0;

	for (uint32_t i = start_pos; i < start_pos + size; ++i)
		if (i & 1)
			r = (r << 2) + (data[i / 2] & 0x3);
		else
			r = (r << 2) + (data[i / 2] >> 4);

	return r;
}

//**********************************************************************************************************
//
//**********************************************************************************************************
void CMappingCore::process_reads_init_stage(reads_bin_t bin)
{
	int64_t dir_exactMatchCounter = 0;
	int64_t rc_exactMatchCounter = 0;
	int64_t dir_and_rc_exactMatchCounter = 0;
	int64_t dir_or_rc_exactMatchCounter = 0;
	int64_t dir_oneMismatchCounter = 0;
	int64_t rc_oneMismatchCounter = 0;

	no_of_Lev_positive = 0;
	no_of_CF_accepted = 0;
	no_of_CF_discarded = 0;
	no_of_CF_accepted_in_malicious_group = 0;
	no_of_CF_discarded_in_malicious_group = 0;

	cmp_lut_ptr = &cmp_lut[0];

	uint32_t bin_id = 0;
	read_id_t id;
	uint32_t size;
	uchar_t *data;
	uchar_t *data_sft;
	uchar_t *old_data = nullptr;
	uint32_t old_data_size = 0;
	uchar_t *rc_data = new uchar_t[max_read_len];
	uchar_t *rc_data_sft = new uchar_t[max_read_len + 1];
	uchar_t *old_matching_part = new uchar_t[max_read_len];
	old_matching_part[0] = 0xff;
	bool newly_read_data = false;
	bool matching_part_is_new = false;
	bool basic_data_part_is_new;

	uint32_t best_error;
	uint32_t prev_read_prefix = ~((uint32_t)0);
	uint32_t *sa_dir_part = nullptr;		// just to make the following code simpler (always delete this array, no ifs necessary)
	uint32_t *sa_rc_part = nullptr;
	uint32_t sa_dir_size, sa_rc_size;

	uint32_t sa_dir_last_index = 0;					//the last index that holds a suffix with a matching part
	uint32_t sa_dir_start_index = 0;					//the first index that holds a suffix with a matching part
	uint32_t sa_rc_last_index = 0;					//the last index that holds a suffix with a matching part
	uint32_t sa_rc_start_index = 0;					//the first index that holds a suffix with a matching part

	uint32_t sa_dir_last_index_for_exact_match = 0;					//the last index that holds a suffix with a basic part
	uint32_t sa_dir_start_index_for_exact_match = 0;					//the first index that holds a suffix with a basic part
	uint32_t sa_rc_last_index_for_exact_match = 0;					//the last index that holds a suffix with a basic part
	uint32_t sa_rc_start_index_for_exact_match = 0;					//the first index that holds a suffix with a basic part

	bool dirExact_match_found = false;
	bool rcExact_match_found = false;

	uint32_t dir_mismatch = 255;
	uint32_t rc_mismatch = 255;

	ref_ptr = reference->GetData();
	ref_size = sa_dir->GetRefSize();

	int64_t no_of_reads = 0;

	sum_of_localizations_dir.clear();
	sum_of_localizations_rc.clear();
	sum_of_localizations_dir.shrink_to_fit();
	sum_of_localizations_rc.shrink_to_fit();

	list_merger_tmp_dir.clear();
	list_merger_tmp_rc.clear();
	list_merger_tmp_dir.shrink_to_fit();
	list_merger_tmp_rc.shrink_to_fit();

	identical_sequences_dir.clear();
	identical_sequences_rc.clear();
	identical_sequences_dir.shrink_to_fit();
	identical_sequences_rc.shrink_to_fit();

	while (reads_deliverer->Pop(id, data, data_sft, size, best_error))//id number was assigned internally
																	  //data_sft - same data as in "data", but shifted by 4 bits
																	  //size - a size of a read in characters
	{
		// Remove too short reads
		if (size < min_read_len)
			continue;

		if (compare_str_to_str(data, min_read_len, old_data))
		{
			basic_data_part_is_new = false;	//along with min_read_len data and old data are the same, 
											//but the rest should be checked

			if (size == old_data_size)
			{
				if (size == min_read_len)
					newly_read_data = false;    //we have another copy of a read of a size of min_read_led
				else
				{
					if (compare_str_to_str(data, size, old_data))
						newly_read_data = false;	//the rest of a read is the same as of previous read
					else
						newly_read_data = true;		//the rest of a read is different
				}
			}
			else
				newly_read_data = true;
		}
		else
		{
			basic_data_part_is_new = true; //along with min_read_len data and old data are different
			newly_read_data = true;
		}
		
		if (newly_read_data)
		{
			//read is new
			// Compute bin id for next substage
			bin_id = get_next_bin_id(data);

			dirExact_match_found = false;
			rcExact_match_found = false;
			old_data = data;
			old_data_size = size;
			dir_match_positions.clear();
			rc_match_positions.clear();

			//			collected_mappings.clear();
			collected_mappings.Clear(max_no_mappings);

			// Compare (prefix_len+sa_prefix_overhead)-symbols long prefix of a read with a previous read
			uint32_t cur_read_prefix = get_read_prefix(data, prefix_len + sa_prefix_overhead);
			//prefix_len - the number of characters that implies a division into bins

			if (cur_read_prefix != prev_read_prefix)
			{
				// Read SA part if necessary
				if (verbosity_level > 1)
					cerr << "Loading SA: " << hex << cur_read_prefix << dec << "   (prefix_len: " << prefix_len << ")\n";
				sa_dir->ReleaseSAPart(sa_dir_part);
				sa_rc->ReleaseSAPart(sa_rc_part);

				sa_dir->GetSAPart(cur_read_prefix, prefix_len + sa_prefix_overhead, sa_dir_part, sa_dir_size);
//				sa_dir_cur_index = 0;
				sa_dir_start_index = 0;
//				start_index = 0;
				sa_dir_last_index = 0;
//				last_index = 0;
				sa_dir_start_index_for_exact_match = 0;
				sa_dir_last_index_for_exact_match = 0;

				sa_rc->GetSAPart(cur_read_prefix, prefix_len + sa_prefix_overhead, sa_rc_part, sa_rc_size);
//				sa_rc_cur_index = 0;
				sa_rc_start_index = 0;
				sa_rc_last_index = 0;
				sa_rc_start_index_for_exact_match = 0;
				sa_rc_last_index_for_exact_match = 0;

#ifdef USE_128b_SSE_CF_AND_DP
				CF_vector_128b_SSE_dir_gen.clear();
				CF_vector_128b_SSE_dir_gen.resize(sa_dir_size, 0);

				CF_vector_128b_SSE_rc_gen.clear();
				CF_vector_128b_SSE_rc_gen.resize(sa_rc_size, 0);
#endif
				prev_read_prefix = cur_read_prefix;
				if (verbosity_level > 2)
					cerr << "Loaded SA : " << hex << cur_read_prefix << dec << "   dir_size: " << sa_dir_size << "   rc_size: " << sa_rc_size << "\n";
			} //end_if(cur_read_prefix != prev_read_prefix)

			reverse_read(data, size, &rc_data, &rc_data_sft, rc_lut);

			//************************ Mapping *****************************

			if (mapping_mode == mapping_mode_t::first)
			{
				dirExact_match_found = findExactMatch_diff_len_version(data, data_sft, size, basic_data_part_is_new, &sa_dir_last_index_for_exact_match, &sa_dir_start_index_for_exact_match, sa_dir_part, sa_dir_size, true);
				rcExact_match_found = findExactMatch_diff_len_version(rc_data, rc_data_sft, size, basic_data_part_is_new, &sa_rc_last_index_for_exact_match, &sa_rc_start_index_for_exact_match, sa_rc_part, sa_rc_size, false);
			}
		}   //end_if(newly_read_data)

			//************************ After mapping ***************************
		if (dirExact_match_found || rcExact_match_found)
		{
			best_error = 0;
			dir_or_rc_exactMatchCounter++;

			if (dirExact_match_found)
			{
				dir_exactMatchCounter++;
				for (uint32_t i = 0; i < dir_match_positions.size(); i++)
					//					results_collector->Push(id, dir_match_positions[i].first, genome_t::direct, 0);
									//					collected_mappings.emplace_back(make_tuple(id, dir_match_positions[i].first, genome_t::direct, 0));
					collected_mappings.Push(dir_match_positions[i].first, genome_t::direct, 0);
			}
			if (rcExact_match_found)
			{
				rc_exactMatchCounter++;
				for (uint32_t i = 0; i < rc_match_positions.size(); i++)
//					results_collector->Push(id, rc_match_positions[i].first, genome_t::rev_comp, 0);
				//					collected_mappings.emplace_back(make_tuple(id, rc_match_positions[i].first, genome_t::rev_comp, 0));
					collected_mappings.Push(rc_match_positions[i].first, genome_t::rev_comp, 0);
			}

			if (dirExact_match_found && rcExact_match_found)
				dir_and_rc_exactMatchCounter++;

			/*			if(collected_mappings.size() > max_no_mappings)
			random_shuffle(collected_mappings.begin(), collected_mappings.end());
			for (auto &x : collected_mappings)
			results_collector->Push(get<0>(x), get<1>(x), get<2>(x), get<3>(x));
			collected_mappings.clear();*/

			while (!collected_mappings.Empty())
			{
				uint64_t x = collected_mappings.PopUnsorted();
				results_collector->Push(id, collected_mappings.DecodePos(x), collected_mappings.DecodeDir(x), collected_mappings.DecodeNoErrors(x));
			}
		}//end_if(dirExact_match_found || rcExact_match_found)

		 //******************************************************************************************************
		 // Neither dirExact_match_found nor rcExact_match_found. Try with one mismatch placed on the second half
		 //******************************************************************************************************
		else
		{
			if (newly_read_data)
			{
				if (!compare_str_to_str(data, predef_match_length, old_matching_part))	//mamy pierwsze wystąpienie matching_part
				{
					string_copy(old_matching_part, data, predef_match_length);
					matching_part_is_new = true;
				}
				else
					matching_part_is_new = false;

				//************************************************************
				// Mapping
				//************************************************************

				dir_mismatch = findMismatches_dir_and_rc(data, data_sft, size, matching_part_is_new, sa_dir_part, sa_dir_size, &sa_dir_last_index, &sa_dir_start_index, 1, &dir_inside_malicious_group, &dir_aux_struct_valid, stage_major);
				rc_mismatch = findMismatches_dir_and_rc(rc_data, rc_data_sft, size, matching_part_is_new, sa_rc_part, sa_rc_size, &sa_rc_last_index, &sa_rc_start_index, 0, &rc_inside_malicious_group, &rc_aux_struct_valid, stage_major);
			}

			//************************************************************
			// After mapping
			//************************************************************
			if ((dir_mismatch <= 1) || (rc_mismatch <= 1))
			{
				if (dir_mismatch <= 1)
				{
					dir_oneMismatchCounter++;
					for (uint32_t i = 0; i < dir_match_positions.size(); i++)
						//						results_collector->Push(id, dir_match_positions[i].first, genome_t::direct, stage_major, 0);
						//						results_collector->Push(id, dir_match_positions[i].first, genome_t::direct, stage_major);
//						results_collector->Push(id, dir_match_positions[i].first, genome_t::direct, dir_mismatch);
					//						collected_mappings.emplace_back(make_tuple(id, dir_match_positions[i].first, genome_t::direct, dir_mismatch));
						collected_mappings.Push(dir_match_positions[i].first, genome_t::direct, dir_mismatch);
				}

				if (rc_mismatch <= 1)
				{
					rc_oneMismatchCounter++;
					for (uint32_t i = 0; i < rc_match_positions.size(); i++)
						//						results_collector->Push(id, rc_match_positions[i].first, genome_t::rev_comp, stage_major, 0);
						//						results_collector->Push(id, rc_match_positions[i].first, genome_t::rev_comp, stage_major);
//						results_collector->Push(id, rc_match_positions[i].first, genome_t::rev_comp, rc_mismatch);
					//						collected_mappings.emplace_back(make_tuple(id, rc_match_positions[i].first, genome_t::rev_comp, rc_mismatch));
						collected_mappings.Push(rc_match_positions[i].first, genome_t::rev_comp, rc_mismatch);
				}
				best_error = MIN(dir_mismatch, rc_mismatch);

				/*				if (collected_mappings.size() > max_no_mappings)
				random_shuffle(collected_mappings.begin(), collected_mappings.end());
				for (auto &x : collected_mappings)
				results_collector->Push(get<0>(x), get<1>(x), get<2>(x), get<3>(x));
				collected_mappings.clear();*/
				while (!collected_mappings.Empty())
				{
					uint64_t x = collected_mappings.PopUnsorted();
					results_collector->Push(id, collected_mappings.DecodePos(x), collected_mappings.DecodeDir(x), collected_mappings.DecodeNoErrors(x));
				}
			}
		} //end else - nither dirExact_match_found nor rcExact_match_found

  	    // Push read for next substage
		bool pass_read = false;

		if (mapping_mode == mapping_mode_t::first)
		{
			if (best_error > 0)
				pass_read = true;
		}
		else if (mapping_mode == mapping_mode_t::second)
		{
			pass_read = true;
		}

		if (pass_read)
		{
			reads_collector->Push(bin_id, id, data, size, best_error);
			no_of_reads++;
		}
	}

	running_stats->AddTotals(STAT_DIR_EXACT_MATCH, dir_exactMatchCounter);
	running_stats->AddTotals(STAT_RC_EXACT_MATCH, rc_exactMatchCounter);
	running_stats->AddTotals(STAT_DIR_AND_RC_EXACT_MATCH, dir_and_rc_exactMatchCounter);
	running_stats->AddTotals(STAT_DIR_OR_RC_EXACT_MATCH, dir_or_rc_exactMatchCounter);

	running_stats->AddTotals(STAT_MISMATCHES_TO_DIR_BASE + stage_id, dir_oneMismatchCounter);
	running_stats->AddTotals(STAT_MISMATCHES_TO_RC_BASE + stage_id, rc_oneMismatchCounter);

	running_stats->AddTotals(STAT_SELECTED_READS_BASE + stage_id, no_of_reads);

	running_stats->AddTotals(STAT_CF_ACCEPTED + stage_id, no_of_CF_accepted);
	running_stats->AddTotals(STAT_CF_DISCARDED + stage_id, no_of_CF_discarded);
	running_stats->AddTotals(STAT_CF_ACCEPTED_IN_MALICIOUS_GROUPS + stage_id, no_of_CF_accepted_in_malicious_group);
	running_stats->AddTotals(STAT_CF_DISCARDED_IN_MALICIOUS_GROUPS + stage_id, no_of_CF_discarded_in_malicious_group);

	running_stats->AddTotals(STAT_LEV_POSITIVE + stage_id, no_of_Lev_positive);

	sa_dir->ReleaseSAPart(sa_dir_part);
	sa_rc->ReleaseSAPart(sa_rc_part);

	delete[] rc_data;
	delete[] rc_data_sft;
	delete[] old_matching_part;
}

//**********************************************************************************************************
// Copy the reads with Ns in prefix
//**********************************************************************************************************
void CMappingCore::copy_reads(reads_bin_t bin)
{
	read_id_t id;
	uchar_t *data;
	uchar_t *data_sft;
	uint32_t size;
	uint32_t best_error;
	int64_t no_of_reads = 0;

	while (reads_deliverer->Pop(id, data, data_sft, size, best_error))
	{
		if (size < min_read_len)		// Remove too short reads
			continue;

		// Compute bid id for next substage
		uint32_t bin_id = get_next_bin_id(data);

		// Push read for next substage
		bool pass_read = false;

		if (stage_minor < stage_major)
			pass_read = true;
		else if (mapping_mode == mapping_mode_t::first)
		{
			if (best_error > stage_major)
				pass_read = true;
		}
		else if (mapping_mode == mapping_mode_t::second)
		{
			if (best_error >= stage_major)
				pass_read = true;
		}

		if (pass_read)
		{
			reads_collector->Push(bin_id, id, data, size, best_error, true);
			no_of_reads++;
		}

		++n_reads;
	}

	running_stats->AddTotals(STAT_SELECTED_READS_BASE + stage_id, no_of_reads);
}

//**********************************************************************************************************
void CMappingCore::complete()
{
	reads_collector->Complete();
	results_collector->Complete();
}

//**********************************************************************************************************
void CMappingCore::process_reads(reads_bin_t bin)
{
	int64_t dir_acceptedMismatchCounter = 0;
	int64_t rc_acceptedMismatchCounter = 0;

	cmp_lut_ptr = &cmp_lut[0];

	read_id_t id;
	uchar_t *data;
	uchar_t *data_sft;
	uchar_t *old_data = nullptr;
	uint32_t old_data_size = 0;
	uchar_t *old_matching_part = new uchar_t[max_read_len];
	old_matching_part[0] = 0xff;

	uint32_t size;
	uint32_t best_error;
	uint32_t prev_read_prefix = ~((uint32_t)0);
	uint32_t *sa_dir_part = nullptr;		// just to make the following code simpler (always delete this array, no ifs necessary)
	uint32_t *sa_rc_part = nullptr;
	uint32_t sa_dir_size, sa_rc_size;

	uint32_t sa_dir_last_index = 0;					// the last index of suffixes' group with an appropriate matching part
	uint32_t sa_dir_start_index = 0;				// the first index of suffixes' group with an appropriate matching part
	uint32_t sa_rc_last_index = 0;					// the last index of suffixes' group with an appropriate matching part
	uint32_t sa_rc_start_index = 0;					// the first index of suffixes' group with an appropriate matching part

	uint32_t bin_id;

	uint32_t dir_mismatch = error_unknown;
	uint32_t rc_mismatch = error_unknown;

	ref_ptr = reference->GetData();
	ref_size = sa_dir->GetRefSize();

	uchar_t *rc_data = new uchar_t[max_read_len];
	uchar_t *rc_data_sft = new uchar_t[max_read_len + 1];
	bool newly_read_data = false;
	bool matching_part_is_new = false;

	uint32_t max_mismatches = 255;
	int64_t no_of_reads = 0;

	sum_of_localizations_dir.clear();
	sum_of_localizations_rc.clear();
	sum_of_localizations_dir.shrink_to_fit();
	sum_of_localizations_rc.shrink_to_fit();

	identical_sequences_dir.clear();
	identical_sequences_rc.clear();
	identical_sequences_dir.shrink_to_fit();
	identical_sequences_rc.shrink_to_fit();

	no_of_CF_accepted = 0;
	no_of_CF_discarded = 0;
	no_of_CF_accepted_in_malicious_group = 0;
	no_of_CF_discarded_in_malicious_group = 0;
	no_of_Lev_positive = 0;

	while (reads_deliverer->Pop(id, data, data_sft, size, best_error))//id number was assigned internally
																	  //data_sft - same data as in "data", but shifted by 4 bits
																	  //size - a size of a read in characters
	{
		if (size == old_data_size)
		{
			if (compare_str_to_str(data, size, old_data))
				newly_read_data = false;
			else
				newly_read_data = true;
		}
		else
			newly_read_data = true;

		if (newly_read_data)
		{
			// read is new
			// Compute bin id for next substage
			bin_id = get_next_bin_id(data);

			old_data = data;
			old_data_size = size;
			dir_match_positions.clear();
			rc_match_positions.clear();

			collected_mappings.Clear(max_no_mappings);

			// Compare (prefix_len+sa_prefix_overhead)-symbols long prefix of read with previous read
			uint32_t cur_read_prefix = get_read_prefix(data, prefix_len + sa_prefix_overhead);
			//prefix_len - the number of characters that implies a division into bins

			if (cur_read_prefix != prev_read_prefix)
			{
				// Read SA part if necessary
				if (verbosity_level > 1)
					cerr << "Loading SA: " << hex << cur_read_prefix << dec << "   (prefix_len: " << prefix_len << ")\n";

				sa_dir->ReleaseSAPart(sa_dir_part);
				sa_rc->ReleaseSAPart(sa_rc_part);

				sa_dir->GetSAPart(cur_read_prefix, prefix_len + sa_prefix_overhead, sa_dir_part, sa_dir_size);
//				sa_dir_cur_index = 0;
				sa_dir_start_index = 0;
				sa_dir_last_index = 0;

				sa_rc->GetSAPart(cur_read_prefix, prefix_len + sa_prefix_overhead, sa_rc_part, sa_rc_size);
//				sa_rc_cur_index = 0;
				sa_rc_start_index = 0;
				sa_rc_last_index = 0;

#ifdef USE_128b_SSE_CF_AND_DP
				CF_vector_128b_SSE_dir_gen.clear();
				CF_vector_128b_SSE_dir_gen.resize(sa_dir_size, 0);

				CF_vector_128b_SSE_rc_gen.clear();
				CF_vector_128b_SSE_rc_gen.resize(sa_rc_size, 0);
#endif
				prev_read_prefix = cur_read_prefix;
				if (verbosity_level > 2)
					cerr << "Loaded SA : " << hex << cur_read_prefix << dec << "   dir_size: " << sa_dir_size << "   rc_size: " << sa_rc_size << "\n";

			} //end_if(cur_read_prefix != prev_read_prefix)

			reverse_read(data, size, &rc_data, &rc_data_sft, rc_lut);

			//-----------------------------------------------------------------------------------------------
			//------- check if "matching part" of just readed data is new -----------------------------------
			//-----------------------------------------------------------------------------------------------
			switch (CUR_OFFSET & 1)
			{
			case 0:			//even CUR_OFFSET
				if (!compare_str_to_str(data + CUR_OFFSET / 2, predef_match_length, old_matching_part))	//old_maching_part does not suit just readed data
				{
					string_copy(old_matching_part, data + CUR_OFFSET / 2, predef_match_length);
					matching_part_is_new = true;
				}
				else
					matching_part_is_new = false;
				break;
			case 1:			//odd CUR_OFFSET
				if (!compare_str_to_str(data_sft + CUR_OFFSET / 2 + 1, predef_match_length, old_matching_part))	//old_maching_part does not suit just readed data
				{
					string_copy(old_matching_part, data_sft + CUR_OFFSET / 2 + 1, predef_match_length);
					matching_part_is_new = true;
				}
				else
					matching_part_is_new = false;
				break;
			} //end_of_switch

			  //-----------------------------------------------------------------------------------------------
			  // map R to dirGenome and rcGenome with at most min(stage_major, best_error) mismatches
			  //-----------------------------------------------------------------------------------------------

			if (mapping_mode == mapping_mode_t::first)
				max_mismatches = (uint32_t)MIN(stage_major * (sensitive_mode ? sensitivity_factor : 1.0), best_error);
			if (mapping_mode == mapping_mode_t::second)
				max_mismatches = (uint32_t)MIN(stage_major * (sensitive_mode ? sensitivity_factor : 1.0), best_error + 1);
			else if (mapping_mode == mapping_mode_t::all)
				max_mismatches = (uint32_t)(stage_major * (sensitive_mode ? sensitivity_factor : 1.0));

			dir_mismatch = findMismatches_dir_and_rc(data, data_sft, size, matching_part_is_new, sa_dir_part, sa_dir_size, &sa_dir_last_index, &sa_dir_start_index, 1, &dir_inside_malicious_group, &dir_aux_struct_valid, max_mismatches);
			rc_mismatch = findMismatches_dir_and_rc(rc_data, rc_data_sft, size, matching_part_is_new, sa_rc_part, sa_rc_size, &sa_rc_last_index, &sa_rc_start_index, 0, &rc_inside_malicious_group, &rc_aux_struct_valid, max_mismatches);

		} //end_if(newly_read_data)
		  //---------------------------------------------------------------------------------------------

		if ((dir_mismatch <= max_mismatches) || (rc_mismatch <= max_mismatches))
		{
			if (dir_mismatch <= max_mismatches)
			{
				dir_acceptedMismatchCounter++;
				for (uint32_t i = 0; i < dir_match_positions.size(); i++)
//					results_collector->Push(id, dir_match_positions[i].first, genome_t::direct, dir_match_positions[i].second);
				//					collected_mappings.emplace_back(make_tuple(id, dir_match_positions[i].first, genome_t::direct, dir_match_positions[i].second));
					collected_mappings.Push(dir_match_positions[i].first, genome_t::direct, dir_match_positions[i].second);
			}

			if (rc_mismatch <= max_mismatches)
			{
				rc_acceptedMismatchCounter++;
				for (uint32_t i = 0; i < rc_match_positions.size(); i++)
//					results_collector->Push(id, rc_match_positions[i].first, genome_t::rev_comp, rc_match_positions[i].second);
				//						collected_mappings.emplace_back(make_tuple(id, rc_match_positions[i].first, genome_t::rev_comp, rc_match_positions[i].second));
					collected_mappings.Push(rc_match_positions[i].first, genome_t::rev_comp, rc_match_positions[i].second);
			}

			/*			 if (collected_mappings.size() > max_no_mappings)
			random_shuffle(collected_mappings.begin(), collected_mappings.end());
			for (auto &x : collected_mappings)
			results_collector->Push(get<0>(x), get<1>(x), get<2>(x), get<3>(x));
			collected_mappings.clear();*/

			while (!collected_mappings.Empty())
			{
				uint64_t x = collected_mappings.PopUnsorted();
				results_collector->Push(id, collected_mappings.DecodePos(x), collected_mappings.DecodeDir(x), collected_mappings.DecodeNoErrors(x));
			}

			if (mapping_mode == mapping_mode_t::first)
			{
				uint32_t tmp_best_error = MIN(dir_mismatch, rc_mismatch);
				best_error = MIN(best_error, tmp_best_error);
			}
			else if (mapping_mode == mapping_mode_t::second)
			{
				uint32_t tmp_best_error = MIN(dir_mismatch, rc_mismatch);
				best_error = MIN(best_error, tmp_best_error);
			}
			else if (mapping_mode == mapping_mode_t::all)
			{
				uint32_t tmp_best_error = MIN(dir_mismatch, rc_mismatch);
				best_error = MIN(best_error, tmp_best_error);
			}
		}  //end if((dir_mismatch <= stage_major) || (rc_mismatch <= stage_major))

		bool pass_read = false;

		if (stage_minor < stage_major)
			pass_read = true;
		else if (mapping_mode == mapping_mode_t::first)
		{
			if (best_error > stage_major)
				pass_read = true;
		}
		else if (mapping_mode == mapping_mode_t::second)
		{
			if (best_error >= stage_major)
				pass_read = true;
		}

		if (pass_read)
		{
			reads_collector->Push(bin_id, id, data, size, best_error);
			no_of_reads++;
		}
	}	//end while

	running_stats->AddTotals(STAT_MISMATCHES_TO_DIR_BASE + stage_id, dir_acceptedMismatchCounter);
	running_stats->AddTotals(STAT_MISMATCHES_TO_RC_BASE + stage_id, rc_acceptedMismatchCounter);
	running_stats->AddTotals(STAT_SELECTED_READS_BASE + stage_id, no_of_reads);

	running_stats->AddTotals(STAT_CF_ACCEPTED + stage_id, no_of_CF_accepted);
	running_stats->AddTotals(STAT_CF_DISCARDED + stage_id, no_of_CF_discarded);
	running_stats->AddTotals(STAT_CF_ACCEPTED_IN_MALICIOUS_GROUPS + stage_id, no_of_CF_accepted_in_malicious_group);
	running_stats->AddTotals(STAT_CF_DISCARDED_IN_MALICIOUS_GROUPS + stage_id, no_of_CF_discarded_in_malicious_group);

	running_stats->AddTotals(STAT_LEV_POSITIVE + stage_id, no_of_Lev_positive);
	//	running_stats->AddTotals(STAT_LEV_NEGATIVE+stage_id, no_of_Lev_negative);

	sa_dir->ReleaseSAPart(sa_dir_part);
	sa_rc->ReleaseSAPart(sa_rc_part);
	delete[] rc_data;
	delete[] rc_data_sft;
	delete[] old_matching_part;
}

//**********************************************************************************************************
void CMappingCore::reverse_read(uchar_t *pattern_ptr, uint32_t pattern_len, uchar_t **rc_pattern_ptr, uchar_t **rc_pattern_ptr_sft, uchar_t *reverse_lut)
{
	uchar_t reverse_lut_for_last_byte[5] = { 0x03, 0x02, 0x01, 0x00, 0x04 };
	uint32_t i, j;

	if (pattern_len & 1)		//the read length is odd
	{
		i = 0; j = pattern_len / 2;
		(*rc_pattern_ptr_sft)[i] = reverse_lut_for_last_byte[(pattern_ptr[j]) >> 4];	// the last half-byte became the first half-byte
		i++; j--;

		for (; i <= pattern_len / 2; i++, j--)
			(*rc_pattern_ptr_sft)[i] = reverse_lut[pattern_ptr[j]];

		//now shift the "shifted data" to the left to get "normal" data
		for (i = 0; i < pattern_len / 2; ++i)
			(*rc_pattern_ptr)[i] = (((*rc_pattern_ptr_sft)[i] & 0xf) << 4) + ((*rc_pattern_ptr_sft)[i + 1] >> 4);

		(*rc_pattern_ptr)[i] = (((*rc_pattern_ptr_sft)[i] & 0xf) << 4);
	}
	else				//the read length is even
	{
		for (i = 0, j = pattern_len / 2 - 1; i < pattern_len / 2; i++, j--)
			(*rc_pattern_ptr)[i] = reverse_lut[pattern_ptr[j]];

		// Shift packed sequence by half of a byte
		(*rc_pattern_ptr_sft)[0] = (*rc_pattern_ptr)[0] >> 4;

		for (uint32_t i = 1; i < pattern_len / 2; ++i)
			(*rc_pattern_ptr_sft)[i] = (((*rc_pattern_ptr)[i - 1] & 0xf) << 4) + ((*rc_pattern_ptr)[i] >> 4);
		(*rc_pattern_ptr_sft)[pattern_len / 2] = ((*rc_pattern_ptr)[pattern_len / 2 - 1] & 0xf) << 4;
	}
}

// EOF