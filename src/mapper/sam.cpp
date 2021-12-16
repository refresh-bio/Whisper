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
#include "sam.h"
#include <algorithm>
#include <utility>
#include "../libs/asmlib.h"
#include "../common/timer.h"
#include "../common/utils.h"
#include "vector_utils.h"

#include <string.h>

#include "LevMyers128.h"
#include "LevMyers256.h"

#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <thread>

#include <string>

//#define LOG_UNMAPPED

// ************************************************************************************
CSamGenerator::CSamGenerator(CParams *params, CObjects *objects, CRefSeqDesc *_ref_seq_desc) :
	params(nullptr),
	reference(nullptr),
#ifdef ENABLE_VCF_VARIANTS
	variant_db(nullptr),
#endif
	f_out{ nullptr, nullptr },
	mem_monitor(nullptr),
	q_map_reads(nullptr),
	mp_fastq_blocks(nullptr),
	mp_ext_cigar(nullptr),
	mp_cigar(nullptr),
	mp_cigar_bin(nullptr),
	mp_mdz(nullptr),
	q_sam_blocks(nullptr),
	ptr_pool(nullptr),
	running_stats(nullptr),
	reads_reader(nullptr),
	reads_readers{ nullptr, nullptr },
	ref_seq_desc(nullptr),
	levMyers64(nullptr),
	levMyers128(nullptr),
	levMyers256(nullptr),
	softClipping(nullptr),
	indelMatching(nullptr),
#ifdef ENABLE_VCF_VARIANTS
	variantMatch(nullptr),
#endif
	insertSizeModel(nullptr),
	progress(nullptr),
	verbosity_level(0),
	id_bits_subgroup(0),
	id_bits_local(0),
	mapping_mode(mapping_mode_t::first),
	mapping_orientation(mapping_orientation_t::forward_reverse),
	max_no_mappings(0),
	mapping_counter_size(0),
	id_bytes(0),
	gzipped_SAM_level(3),
	store_BAM(false),
	data_srt(nullptr),
	data_srt_end(nullptr),
	max_mate_distance(0),
	max_mate_edit_distance_frac(0),
	tmp_read_sequence(nullptr),
	tmp_ref_sequence(nullptr),
	sam_part(nullptr),
	bam_part(nullptr),
	stored_mapped(0),
	stored_unmapped(0),
	read_bytes(0),
	stat_mapped_se_single(0),
	stat_mapped_se_none(0),
	stat_mapped_pe_pair_unq(0),
	stat_mapped_pe_pair_more(0),
	stat_mapped_pe_independent_both(0),
	stat_mapped_pe_independent_single(0),
	stat_mapped_pe_independent_none(0),
	stat_mapped_pe_independent_errors(0),
	stat_mapped_pe_errors(0),
	stat_refined_short_indels_lev(0),
	stat_refined_short_indels_mapping(0)
{
	instruction_set_t instruction_set;
	int x = instrset_detect();
	if (x >= 0 && x <= 8)
		instruction_set = (instruction_set_t)x;
	else if (x < 0)
		instruction_set = instruction_set_t::none;
	else
		instruction_set = instruction_set_t::avx2;
	
	reference = objects->reference;
#ifdef ENABLE_VCF_VARIANTS
	variant_db = objects->variant_db;
#endif

	mem_monitor = objects->mem_monitor;

	q_map_reads     = objects->q_map_reads;
	mp_fastq_blocks = objects->mp_fastq_blocks;
	q_sam_blocks    = objects->q_sam_blocks;
	mp_ext_cigar    = objects->mp_ext_cigar;
	mp_cigar		= objects->mp_cigar;
	mp_cigar_bin    = objects->mp_cigar_bin;
	mp_mdz			= objects->mp_mdz;
	ptr_pool = objects->ptr_pool;

	progress = objects->progress;

	running_stats = objects->running_stats;

	this->params = params;

	fill(stat_read_len.begin(), stat_read_len.end(), 0u);
	fill(stat_read_len_wo_Ns.begin(), stat_read_len_wo_Ns.end(), 0u);
	fill(stat_read_Ns.begin(), stat_read_Ns.end(), 0u);
	
	verbosity_level      = params->verbosity_level;
	id_bits_subgroup     = params->id_bits_subgroup;
	id_bits_local        = params->id_bits_local;

	mapping_mode		 = params->mapping_mode;
	mapping_orientation  = params->mapping_orientation;
	max_no_mappings		 = params->max_no_mappings;
	mapping_counter_size = params->mapping_counter_size;

	max_mate_distance			= params->max_mate_distance;
	max_mate_edit_distance_frac = params->max_mate_edit_distance_frac;

	id_bytes             = BITS2BYTES(params->id_bits_total);

	gzipped_SAM_level = params->gzipped_SAM_level;
	store_BAM = params->store_BAM;

	if (store_BAM)
	{
		sam_part = nullptr;
		bam_part = new CBamPart(params, objects);
	}
	else
	{
		sam_part = new CSamPart(params, objects);
		bam_part = nullptr;
	}

	ref_seq_desc         = _ref_seq_desc;

	reads_reader         = new CReadsReader(params, objects, params->is_fasta);
	reads_readers[0]     = new CReadsReader(params, objects, params->is_fasta);
	reads_readers[1]     = new CReadsReader(params, objects, params->is_fasta);

	stored_mapped   = 0;
	
	read_bytes      = 0;

	scoring_t scoring(
		params->match_score,
		params->mismatch_score,
		params->mismatch_score,
//		params->gap_open,
//		params->gap_extend,
		params->gap_ins_open,
		params->gap_ins_extend,
		params->gap_del_open,
		params->gap_del_extend,
		params->clipping_score);
	
	insertSizeModel = new InsertSizeModel(
		max_mate_distance / 2, max_mate_distance / 8, 
		100, 10000,
		params->high_confidence_sigmas, params->clipping_distance_sigmas, params->penalty_saturation_sigmas, scoring);

	int max_text_len = 16384;

	levMyers64 = new LevMyers64(params->max_read_len, max_text_len, 0);

	switch (instruction_set) {
	case instruction_set_t::none:
	case instruction_set_t::sse:
		throw new std::runtime_error("SSE2 extensions required!");
	case instruction_set_t::sse2:
	case instruction_set_t::sse3:
	case instruction_set_t::sse3s:
	case instruction_set_t::sse41:
	case instruction_set_t::sse42:
		levMyers128 = new LevMyers128<instruction_set_t::sse2>(max_text_len, 0);
		levMyers256 = new LevMyers64(params->max_read_len, max_text_len, 0);
		break;
	case instruction_set_t::avx:
		levMyers128 = new LevMyers128<instruction_set_t::avx>(max_text_len, 0);
		levMyers256 = new LevMyers256<instruction_set_t::avx>(max_text_len, 0);
		break;
	case instruction_set_t::avx2:
		levMyers128 = new LevMyers128<instruction_set_t::avx2>(max_text_len, 0);
		levMyers256 = new LevMyers256<instruction_set_t::avx2>(max_text_len, 0);
		break;
	}

	levMyers64->setReference(reference->GetData(), 0, 0);
	levMyers128->setReference(reference->GetData(), 0, 0);
	levMyers256->setReference(reference->GetData(), 0, 0);
	
	softClipping = new CSoftClipping(params->max_read_len, max_text_len, params->min_approx_indel_len, params->max_approx_indel_len, params->max_approx_indel_mismatches);
	softClipping->setReference(reference->GetData(), 0, 0);

	indelMatching = new CIndelMatching(params);
	indelMatching->SetReference(reference->GetData(), 0, 0);

#ifdef ENABLE_VCF_VARIANTS
	variantMatch = new CVariantMatch();
#endif

	tmp_read_sequence = new uchar_t[params->max_read_len + 1ull];

	tmp_ref_sequence  = new uchar_t[max_text_len + 1ull];

	mapping_desc[0].reserve(max_no_mappings + 1ull);
	mapping_desc[1].reserve(max_no_mappings + 1ull);
}

// ************************************************************************************
CSamGenerator::~CSamGenerator()
{
	delete reads_reader;
	delete reads_readers[0];
	delete reads_readers[1];

	delete levMyers256; 
	delete levMyers128;
	delete levMyers64; 

	delete softClipping;
	delete indelMatching;

#ifdef ENABLE_VCF_VARIANTS
	delete variantMatch;
#endif

	delete insertSizeModel;

	delete[] tmp_read_sequence;
	delete[] tmp_ref_sequence;

	if(sam_part)
		delete sam_part;
	if (bam_part)
		delete bam_part;

#ifdef COLLECT_STATS
	running_stats->AddValues(1000002, read_bytes);
#endif
}

// ************************************************************************************
void CSamGenerator::operator()()
{
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	stat_mapped_se_single    = 0;
	stat_mapped_se_none      = 0;
	stat_mapped_pe_pair_unq  = 0;
	stat_mapped_pe_pair_more = 0;

	stat_mapped_pe_errors = 0;

	stat_mapped_pe_independent_both = 0;
	stat_mapped_pe_independent_single = 0;
	stat_mapped_pe_independent_none = 0;
	stat_mapped_pe_independent_errors = 0;
	
	stat_mapped_pe_histo.resize(50, 0);
	stat_mapped_pe_histo_clipped.resize(50, 0);
	stat_mapped_pe_method.resize(MatchingMethod::size(), 0);
	
	stat_refined_short_indels_lev = 0;
	stat_refined_short_indels_mapping = 0;

#ifdef LOG_UNMAPPED
	uint32_t thread_id = (uint32_t)hash<thread::id>()(this_thread::get_id());

	for (int i = 0; i < 2; ++i)
	{
		string fn = string("unmapped_") + to_string(thread_id) + "_" + to_string(i + 1) + ".fastq";
		f_out[i] = fopen(fn.c_str(), "wb");
		setvbuf(f_out[i], nullptr, _IOFBF, 32 << 20);
	}
#endif

	while(!q_map_reads->IsCompleted())
	{
// Ommit some bins - just for debugging
#if 0
		if (q_map_reads->Pop(mapped_reads))
		{
			int bin_id = (mapped_reads.fastq_blocks.front().id_range >> (id_bits_subgroup + id_bits_local));
			cerr << bin_id << "     " << endl; fflush(stderr);

			if (bin_id == 298)
//			if (bin_id == 1)
			{
				if (verbosity_level >= 2)
					cerr << "SAM generation of results bin no. : " << (mapped_reads.fastq_blocks.front().id_range >> (id_bits_subgroup + id_bits_local)) << "\n";

				try {
					sort_mapping_results();
				}
				catch (...)
				{
					cerr << "\nError in sort\n";
					fflush(stderr);
				}

				try {
					if (mapped_reads.single_end)
						process_group_se();		// Single end reads
					else
						process_group_pe(bin_id);		// Paired end eads
				}
				catch (...)
				{
					cerr << "\nError in processing bin\n";
					fflush(stderr);
				}
			}
			else
			{
				vector<fastq_block_t> fastq_blocks[2];
//				ptrs = nullptr;

				// Split FASTQ blocks into left and right
				for (auto &p : mapped_reads.fastq_blocks)
					fastq_blocks[p.id_range & 1].push_back(p);

				uint32_t no_fastq_blocks = (uint32_t)fastq_blocks[0].size();
				for (uint32_t i = 0; i < no_fastq_blocks; ++i)
				{
					mp_fastq_blocks->Free(fastq_blocks[0][i].data);
					mp_fastq_blocks->Free(fastq_blocks[1][i].data);
				}

				mapped_reads.fastq_blocks.clear();
			}

			progress->Step(1);

			free_reads();
		}
#else
		if(q_map_reads->Pop(mapped_reads))
		{
			if(verbosity_level >= 2)
				cerr << "SAM generation of results bin no. : " << (mapped_reads.fastq_blocks.front().id_range >> (id_bits_subgroup + id_bits_local)) << "\n";

			sort_mapping_results();

			if(mapped_reads.single_end)
				process_group_se();		// Single end reads
			else
				process_group_pe((mapped_reads.fastq_blocks.front().id_range >> (id_bits_subgroup + id_bits_local)));		// Paired end eads

			progress->Step(1);

			free_reads();
		}
#endif

	}

	q_sam_blocks->MarkCompleted();

	thr_watch.StopTimer();

#ifdef COLLECT_STATS
	running_stats->AddValues(STAT_TIME_THR_SAM_GENERATOR, thr_watch.GetElapsedTime());

	running_stats->AddTotals(STAT_MAPPED_SE_SINGLE,    stat_mapped_se_single);
	running_stats->AddTotals(STAT_MAPPED_SE_NONE,      stat_mapped_se_none);
	running_stats->AddTotals(STAT_MAPPED_PE_PAIR_UNQ,  stat_mapped_pe_pair_unq);
	running_stats->AddTotals(STAT_MAPPED_PE_PAIR_MORE, stat_mapped_pe_pair_more);

	running_stats->AddTotals(STAT_MAPPED_PE_ERRORS, stat_mapped_pe_errors);
	
	running_stats->AddTotals(STAT_MAPPED_PE_INDEPENDENT_BOTH, stat_mapped_pe_independent_both);
	running_stats->AddTotals(STAT_MAPPED_PE_INDEPENDENT_SINGLE, stat_mapped_pe_independent_single);
	running_stats->AddTotals(STAT_MAPPED_PE_INDEPENDENT_NONE, stat_mapped_pe_independent_none);
	running_stats->AddTotals(STAT_MAPPED_PE_INDEPENDENT_ERRORS, stat_mapped_pe_independent_errors);

	running_stats->AddTotals(STAT_SHORT_INDEL_REFINEMENTS_LEV, stat_refined_short_indels_lev);
	running_stats->AddTotals(STAT_SHORT_INDEL_REFINEMENTS_MAPPING, stat_refined_short_indels_mapping);

	for (int i = 0; i < (int)stat_mapped_pe_method.size(); ++i) {
		running_stats->AddTotals(STAT_MAPPED_PE_METHOD + i, stat_mapped_pe_method[i]);
	}

	for (int i = 0; i < (int) stat_mapped_pe_histo.size(); ++i) {
		running_stats->AddTotals(STAT_MAPPED_PE_HISTO + i, stat_mapped_pe_histo[i]);
	}

	for (int i = 0; i < (int) stat_mapped_pe_histo_clipped.size(); ++i) {
		running_stats->AddTotals(STAT_MAPPED_PE_HISTO_CLIPPED + i, stat_mapped_pe_histo_clipped[i]);
	}

	for (int i = 0; i <= (int)STAT_MAX_READ_LEN; ++i)
	{
		running_stats->AddTotals(STAT_READS_LEN + i, stat_read_len[i]);
		running_stats->AddTotals(STAT_READS_NS + i, stat_read_Ns[i]);
		running_stats->AddTotals(STAT_READS_LEN_WO_NS + i, stat_read_len_wo_Ns[i]);
	}
#endif

#ifdef LOG_UNMAPPED
	for (int i = 0; i < 2; ++i)
		fclose(f_out[i]);
#endif
}

// ************************************************************************************
void CSamGenerator::sort_mapping_results()
{
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	uint64_t mapped_size = mapped_reads.results.size;
	uint64_t max_no_map_res = mapped_size / candidate_mapping_t::min_size + 8;
	uint64_t last_idx = max_no_map_res;
	uchar_t *data = mapped_reads.results.data;
	read_id_t prev_id = empty_read_id;
	read_id_t id;
	uint64_t no_mappings_in_curr_read = 0;

	if(max_no_map_res > ptrs.capacity())
		ptrs.reserve((size_t) (max_no_map_res * 1.1));		// 10% overhead to reduce the number of reallocations for future bins

	ptrs.resize(max_no_map_res);

	// Localize all read ids
	for(uint64_t i = 0; i < mapped_size;)
	{
		--last_idx;
		ptrs[last_idx].first = i;
		uint64_t no_mappings;
		LoadUInt(data+i+id_bytes, no_mappings, mapping_counter_size);
		uint32_t no_stored_mappings = (uint32_t) no_mappings;
		size_t desc_size = 0;
		uchar_t* ptr = data + i + id_bytes + mapping_counter_size;

		for (uint32_t j = 0; j < no_stored_mappings; ++j)
			desc_size += candidate_mapping_t::check_rec_size(ptr + desc_size);

		ptrs[last_idx].second = desc_size;
		i += id_bytes + mapping_counter_size + desc_size;
	}

	if (verbosity_level >= 2)
		cerr << "Sorting " << max_no_map_res-last_idx << " results\n";
	
	// Sort pointers to read mapping results
	sort(ptrs.data()+last_idx, ptrs.data()+max_no_map_res, [&](pair<uint64_t, uint64_t> x, pair<uint64_t, uint64_t> y){
		for(uint32_t i = 0; i < id_bytes; ++i)
			if(data[x.first+i] < data[y.first+i])
				return true;
			else if(data[x.first+i] > data[y.first+i])
				return false;
		return false;
	});

	// Reorganizing mapped results according to read id
	uint64_t res_idx = 0;
	data_srt = (uchar_t*) ptrs.data();
	
	uint64_t cur_mapping_counter_pos;

	for(uint64_t i = last_idx; i < max_no_map_res; ++i)
	{
		uint64_t tmp;
		uint64_t no_mappings_in_curr_part;

		LoadUInt(data+ptrs[i].first, tmp, id_bytes);
		id = tmp;
		LoadUInt(data+ptrs[i].first+id_bytes, no_mappings_in_curr_part, mapping_counter_size);
			
		if(id != prev_id)
		{
			copy_n(data+ptrs[i].first, id_bytes + mapping_counter_size, data_srt + res_idx);
			res_idx += (uint64_t) id_bytes + mapping_counter_size;
			cur_mapping_counter_pos = res_idx - mapping_counter_size;
			prev_id = id;
			no_mappings_in_curr_read = 0;
		}

		no_mappings_in_curr_read += no_mappings_in_curr_part;
		StoreUInt(data_srt + cur_mapping_counter_pos, no_mappings_in_curr_read, mapping_counter_size);

		copy_n(data+ptrs[i].first + id_bytes + mapping_counter_size, 
			ptrs[i].second,
			data_srt + res_idx);
		res_idx += ptrs[i].second;
	}

	data_srt_end = data_srt + res_idx;

	thr_watch.StopTimer();

#ifdef COLLECT_STATS
	running_stats->AddValues(STAT_TIME_THR_SAM_SORTING, thr_watch.GetElapsedTime());
#endif
}

// ************************************************************************************
void CSamGenerator::push_mapped_part()
{
	if (store_BAM)
	{
		q_sam_blocks->Push(bam_part->GetBlock());
	}
	else
	{
		if (gzipped_SAM_level == 0)
			q_sam_blocks->Push(sam_part->GetBlock());
		else
		{
#ifndef _DEBUG
//			uchar_t *gz_buffer;
			//		mp_sam_parts->Reserve(gz_buffer);
			//		size_t comp_size = gzip->Compress(mapped_part, mapped_part_pos, gz_buffer, mapped_part_size);
			//		q_sam_blocks->Push(sam_block_t(gz_buffer, comp_size, sam_results_t::mapped));
			//		mp_sam_parts->Free(mapped_part);
#endif
		}
	}
}

// ************************************************************************************
// Process group of mapping results for single end reads
void CSamGenerator::process_group_se()
{
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	uchar_t *data_ptr = data_srt;

	for(auto &p : mapped_reads.fastq_blocks)
	{
		uchar_t *id, *sequence, *plus, *quality;
		read_id_t read_id;
		string ref_seq_name;
		uint32_t id_len, sequence_len, plus_len, quality_len;

		if (store_BAM)
			bam_part->Reserve();
		else
			sam_part->Reserve();
	
		read_id = p.id_range;

		read_bytes += p.size;

		reads_reader->SetBlock(p.data, p.size);
		adjust_pos_in_results(read_id, data_ptr);
		while(reads_reader->Pop(id, sequence, plus, quality, id_len, sequence_len, plus_len, quality_len))
		{
			if(find_id(read_id, data_ptr, data_srt_end))
			{
				store_mapped_read(id, sequence, plus, quality, id_len, sequence_len, plus_len, quality_len, data_ptr);
				stat_mapped_se_single++;
			}
			else
			{
				store_unmapped_read(id, sequence, plus, quality, id_len, sequence_len, plus_len, quality_len, data_ptr);
				stat_mapped_se_none++;
			}
			
			read_id += 2;
		}

		push_mapped_part();
		mp_fastq_blocks->Free(p.data);
	}

	mapped_reads.fastq_blocks.clear();

	thr_watch.StopTimer();

#ifdef COLLECT_STATS
	running_stats->AddValues(STAT_TIME_THR_SAM_PROCESSING, thr_watch.GetElapsedTime());
#endif
}

// ************************************************************************************
// Process group of mapping results for paired end reads
void CSamGenerator::process_group_pe(uint64_t group_id)
{
	CThreadWatch thr_watch;
	thr_watch.StartTimer();

#ifdef MEASURE_PP_TIMES
	fill_n(part_times.begin(), 16, 0.0);
#endif

	uchar_t *data_ptr = data_srt;
	vector<fastq_block_t> fastq_blocks[2];

	// Split FASTQ blocks into left and right
	for(auto &p : mapped_reads.fastq_blocks)
		fastq_blocks[p.id_range & 1].push_back(p);

	for(uint32_t i = 0; i < 2; ++i)
		sort(fastq_blocks[i].begin(), fastq_blocks[i].end(), [](fastq_block_t x, fastq_block_t y) {
			return x.id_range < y.id_range;
		});

	uint32_t no_fastq_blocks = (uint32_t) fastq_blocks[0].size();

	insertSizeModel->reset();

	// Process pairs of blocks
	for(uint32_t i = 0; i < no_fastq_blocks; ++i)
	{
		uchar_t *id[2], *sequence[2], *plus[2], *quality[2];
		read_id_t read_id[2];
		string ref_seq_name[2];
		uint32_t id_len[2], sequence_len[2], plus_len[2], quality_len[2];

		if (store_BAM)
			bam_part->Reserve();
		else
			sam_part->Reserve();

		read_id[0] = fastq_blocks[0][i].id_range;
		read_id[1] = fastq_blocks[1][i].id_range;

		read_bytes += fastq_blocks[0][i].size;
		read_bytes += fastq_blocks[1][i].size;

		reads_readers[0]->SetBlock(fastq_blocks[0][i].data, fastq_blocks[0][i].size);
		reads_readers[1]->SetBlock(fastq_blocks[1][i].data, fastq_blocks[1][i].size);

		adjust_pos_in_results(read_id[0], data_ptr);
		uchar_t* ptr_copy = data_ptr;

		// train a model on the first block
		if (i == 0) {
			bool continueTraining = true;
			while (continueTraining) {
				bool data_present = true;
				for (uint32_t j = 0; j < 2; ++j)
					data_present &= reads_readers[j]->Pop(id[j], sequence[j], plus[j], quality[j], id_len[j], sequence_len[j], plus_len[j], quality_len[j]);

				if (!data_present)
					break;

				continueTraining = update_model(read_id[0], id, sequence, id_len, sequence_len, data_ptr);
				
				read_id[0] += 2;
				read_id[1] += 2;
			}

			// restore position
			read_id[0] = fastq_blocks[0][i].id_range;
			read_id[1] = fastq_blocks[1][i].id_range;

			data_ptr = ptr_copy;
			adjust_pos_in_results(read_id[0], data_ptr);
			reads_readers[0]->Restart();
			reads_readers[1]->Restart();
		}

		while(true)
		{			
			bool data_present = true;
			for(uint32_t j = 0; j < 2; ++j)
				data_present &= reads_readers[j]->Pop(id[j], sequence[j], plus[j], quality[j], id_len[j], sequence_len[j], plus_len[j], quality_len[j]);

			if(!data_present)
				break;

			store_mapped_pair_reads(read_id[0], id, sequence, plus, quality, id_len, sequence_len, plus_len, quality_len, data_ptr);

			for (int j = 0; j < 2; ++j)
			{
				stat_read_len[sequence_len[j]]++;
				auto no_Ns = count(sequence[j], sequence[j] + sequence_len[j], 'N');
				stat_read_Ns[no_Ns]++;
				stat_read_len_wo_Ns[sequence_len[j] - no_Ns]++;
			}

			read_id[0] += 2;
			read_id[1] += 2;
		}

		push_mapped_part();
		mp_fastq_blocks->Free(fastq_blocks[0][i].data);
		mp_fastq_blocks->Free(fastq_blocks[1][i].data);
	}
	
	mapped_reads.fastq_blocks.clear();

	thr_watch.StopTimer();

#ifdef COLLECT_STATS
	running_stats->AddValues(STAT_TIME_THR_SAM_PROCESSING, thr_watch.GetElapsedTime());

#ifdef MEASURE_PP_TIMES
	string prefix = "Postprocessing - " + to_string(group_id) + " - ";

	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 0, prefix + "Mappings loading", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 1, prefix + "Mappings sorting", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 2, prefix + "Find pairs: Mapping close", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 3, prefix + "Find pairs: Variants", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 4, prefix + "Find pairs: LevMyers - low range", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 5, prefix + "Find pairs: Clipping - low range", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 6, prefix + "Find paits: LevMyers - high range", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 7, prefix + "Find pairs: Clipping - high range", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 8, prefix + "Find pairs: Mapping distant", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 9, prefix + "Sorting SE", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 10, prefix + "Cigar preparation", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 11, prefix + "Find single", running_stats_t::totals);
	running_stats->Register(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + 12, prefix + "Store results", running_stats_t::totals);

	for (int i = 0; i < 13; ++i)
		running_stats->AddTotals(STAT_TIME_THR_PP_PARTS + (uint32_t) group_id * 32 + i, part_times[i]);
#endif
#endif
}

// ************************************************************************************
void CSamGenerator::free_reads()
{
	ptr_pool->Release(mapped_reads.results.data);
	ptrs.clear();

	mem_monitor->Decrease(mapped_reads.results.size*2);
}

// ************************************************************************************
// Check whether for the current read id (from FASTQ file) there are any mapping results
bool CSamGenerator::find_id(read_id_t read_id, uchar_t *data_ptr, uchar_t *data_srt_end)
{
	uint64_t tmp;
	read_id_t id;

	if(data_ptr >= data_srt_end)
		return false;

	LoadUInt(data_ptr, tmp, id_bytes);
	id = tmp;

	if(read_id == id)
		return true;
	if(read_id < id)
		return false;

	return false;
}

// ************************************************************************************
// Check whether for the current pair (read_id, read_id+1) (from FASTQ file) there are mapping results for any of the reads
bool CSamGenerator::find_id_pair(read_id_t read_id, uchar_t *data_ptr, uchar_t *data_srt_end)
{
	uint64_t tmp;
	read_id_t id;

	if(data_ptr >= data_srt_end)
		return false;

	// Check 1st id
	LoadUInt(data_ptr, tmp, id_bytes);
	id = tmp;

	return (id == read_id || id == read_id+1);

	if(read_id == id)
		return true;

	// Check 2nd id
	LoadUInt(data_ptr, tmp, id_bytes);
	id = tmp;

	if(read_id+1 == id)
		return false;

	return true;
}

// ************************************************************************************
// Extracts CIGAR and MDZ from internal extended cigar format
int CSamGenerator::extract_cigar_and_mdz(const uchar_t * ext_cigar, uchar_t * cigar, uchar_t * mdz)
{
	// auxiliary enum type 
	enum class MappingState {Match, Mismatch, Insertion, Deletion, Begin, End, Clipping};
	
	// iterate over all ext_cigar symbols
	const uchar_t* ext_c = ext_cigar;
	uchar_t* cigar_c = cigar;
	uchar_t* mdz_c = mdz;
	
	if (ext_cigar) {
		int counters[7] = { 0 }; // counters for all states
		int cnt_match_mismatch = 0; // counter for continuous match-mismatch state
		int cnt_match_ignoring_insert = 0; // counter for continuous match-insertions
		MappingState currentState = MappingState::Begin;
		MappingState prevState = MappingState::Begin;

		while (currentState != MappingState::End) {
			
			// detect states and fill MDZ

			if (*ext_c == '$') {
				// clipping 
				currentState = MappingState::Clipping;
				++ext_c;
			} else if (*ext_c == 0) {
				// end state
				currentState = MappingState::End;
				// write number of previous matches
				if (cnt_match_ignoring_insert > 0) {
					mdz_c += CNumericConversions::Int2PChar(cnt_match_ignoring_insert, mdz_c);
				}
			} else if (*ext_c == '.') {
				// match - do nothing
				++ext_c;
				currentState = MappingState::Match;
			}
			else if (*ext_c == '#') {
				// insertion in read - do nothing
				++ext_c;
				++ext_c;
				currentState = MappingState::Insertion;
			}
			else {
				// mismatch or deletion - write number of previous matches
				mdz_c += CNumericConversions::Int2PChar(cnt_match_ignoring_insert, mdz_c);
				cnt_match_ignoring_insert = 0;

				if (*ext_c == '^') {
					// deletion in read - write two elements in MDZ
					*mdz_c++ = *ext_c++;
					*mdz_c++ = *ext_c++;
					currentState = MappingState::Deletion;
				}
				else {
					// put symbol in MDZ
					*mdz_c++ = *ext_c++;
					currentState = MappingState::Mismatch;
				}
			}

			// update counters
			++counters[static_cast<int>(currentState)];
			
			if (currentState == MappingState::Match) 
			{
				++cnt_match_mismatch;
				++cnt_match_ignoring_insert;
			} 
			else if (currentState == MappingState::Mismatch) 
			{
				++cnt_match_mismatch;
			}

			
			// state change (excluding transition from Begin, including transition to End)
			if ((currentState != prevState) && (prevState != MappingState::Begin)) {

				// clipping -> X 
				if (prevState == MappingState::Clipping) {
					cigar_c += CNumericConversions::Int2PChar(counters[static_cast<int>(MappingState::Clipping)], cigar_c);
					*cigar_c++ = 'S';
				}

				// insertion -> X
				else if (prevState == MappingState::Insertion) {
					cigar_c += CNumericConversions::Int2PChar(counters[static_cast<int>(MappingState::Insertion)], cigar_c);
					*cigar_c++ = 'I';
				}
				// deletion -> X
				else if (prevState == MappingState::Deletion) {
					cigar_c += CNumericConversions::Int2PChar(counters[static_cast<int>(MappingState::Deletion)], cigar_c);
					*cigar_c++ = 'D';
				}
				// match -> X
				else if (prevState == MappingState::Match) {
					// match -> not mismatch
					if (currentState != MappingState::Mismatch) {
						cigar_c += CNumericConversions::Int2PChar(cnt_match_mismatch, cigar_c);
						*cigar_c++ = 'M';
						cnt_match_mismatch = 0;
					}
				}
				// mismatch -> not match
				else if (prevState == MappingState::Mismatch) {
					if (currentState != MappingState::Match) {
						cigar_c += CNumericConversions::Int2PChar(cnt_match_mismatch, cigar_c);
						*cigar_c++ = 'M';
						cnt_match_mismatch = 0;
					}
				}

				counters[static_cast<int>(prevState)] = 0;
			}

			prevState = currentState;
		}
	}

	*cigar_c = 0;
	*mdz_c = 0;
	
	return (int) (cigar_c - cigar);
}

// ************************************************************************************
// Extracts binary CIGAR and MDZ from internal extended cigar format
int CSamGenerator::extract_cigarbin_and_mdz(const uchar_t * ext_cigar, uint32_t * cigar_bin, uchar_t * mdz)
{
	// auxiliary enum type 
	enum class MappingState { Match, Mismatch, Insertion, Deletion, Begin, End, Clipping };
	enum class OperationCode { M_op, I_op, D_op, N_op, S_op, H_op, P_op};

	// iterate over all ext_cigar symbols
	const uchar_t* ext_c = ext_cigar;
	uint32_t* cigar_c = cigar_bin;
	uchar_t* mdz_c = mdz;

	if (ext_cigar) {
		int counters[7] = { 0 }; // counters for all states
		int cnt_match_mismatch = 0; // counter for continuous match-mismatch state
		int cnt_match_ignoring_insert = 0; // counter for continuous match-insertions
		MappingState currentState = MappingState::Begin;
		MappingState prevState = MappingState::Begin;

		while (currentState != MappingState::End) {

			// detect states and fill MDZ

			if (*ext_c == '$') {
				// clipping 
				currentState = MappingState::Clipping;
				++ext_c;
			}
			else if (*ext_c == 0) {
				// end state
				currentState = MappingState::End;
				// write number of previous matches
				if (cnt_match_ignoring_insert > 0) {
					mdz_c += CNumericConversions::Int2PChar(cnt_match_ignoring_insert, mdz_c);
				}
			}
			else if (*ext_c == '.') {
				// match - do nothing
				++ext_c;
				currentState = MappingState::Match;
			}
			else if (*ext_c == '#') {
				// insertion in read - do nothing
				++ext_c;
				++ext_c;
				currentState = MappingState::Insertion;
			}
			else {
				// mismatch or deletion - write number of previous matches
				mdz_c += CNumericConversions::Int2PChar(cnt_match_ignoring_insert, mdz_c);
				cnt_match_ignoring_insert = 0;

				if (*ext_c == '^') {
					// deletion in read - write two elements in MDZ
					*mdz_c++ = *ext_c++;
					*mdz_c++ = *ext_c++;
					currentState = MappingState::Deletion;
				}
				else {
					// put symbol in MDZ
					*mdz_c++ = *ext_c++;
					currentState = MappingState::Mismatch;
				}
			}

			// update counters
			++counters[static_cast<int>(currentState)];

			if (currentState == MappingState::Match) {
				++cnt_match_mismatch;
				++cnt_match_ignoring_insert;
			}
			else if (currentState == MappingState::Mismatch) {
				++cnt_match_mismatch;
			}


			// state change (excluding transition from Begin, including transition to End)
			if ((currentState != prevState) && (prevState != MappingState::Begin)) {

				// clipping -> X 
				if (prevState == MappingState::Clipping) {
					*cigar_c++ = counters[static_cast<int>(MappingState::Clipping)] << 4 | static_cast<int>(OperationCode::S_op);
				}

				// insertion -> X
				else if (prevState == MappingState::Insertion) {
					*cigar_c++ = counters[static_cast<int>(MappingState::Insertion)] << 4 | static_cast<int>(OperationCode::I_op);
				}
				// deletion -> X
				else if (prevState == MappingState::Deletion) {
					*cigar_c++ = counters[static_cast<int>(MappingState::Deletion)] << 4 | static_cast<int>(OperationCode::D_op);
				}
				// match -> X
				else if (prevState == MappingState::Match) {
					// match -> not mismatch
					if (currentState != MappingState::Mismatch) {
						*cigar_c++ = cnt_match_mismatch << 4 | static_cast<int>(OperationCode::M_op);
						cnt_match_mismatch = 0;
					}
				}
				// mismatch -> not match
				else if (prevState == MappingState::Mismatch) {
					if (currentState != MappingState::Match) {
						*cigar_c++ = cnt_match_mismatch << 4 | static_cast<int>(OperationCode::M_op);
						cnt_match_mismatch = 0;
					}
				}

				counters[static_cast<int>(prevState)] = 0;
			}

			prevState = currentState;
		}
	}

	*mdz_c = 0;

	return (int) (cigar_c - cigar_bin);
}

// ************************************************************************************
// Adjust position in results
bool CSamGenerator::adjust_pos_in_results(read_id_t read_id, uchar_t *&data_ptr)
{
	data_ptr = data_srt;

	uint64_t tmp;
	uint64_t no_mappings;
	read_id_t id;

	while(data_ptr < data_srt_end)
	{
		LoadUInt(data_ptr, tmp, id_bytes);
		id = tmp;
		LoadUInt(data_ptr+id_bytes, no_mappings, mapping_counter_size);

		if(read_id <= id)
			break;

		data_ptr += (uint64_t) id_bytes + mapping_counter_size;

		uint64_t no_stored_mappings = no_mappings;

		for(uint32_t i = 0; i < no_stored_mappings; ++i)
			data_ptr = candidate_mapping_t::skip(data_ptr);
	}

	return true;
}

// ************************************************************************************
// Computes the length of the first term of id (to the first white space)
inline uint32_t CSamGenerator::trim_id(uchar_t *id, uint32_t id_len)
{
	for(uint32_t r = 0; r < id_len; ++r)
		if(id[r] == ' ' || id[r] == '\t' || id[r] == '\n')
			return r;

	return id_len;
}

// ************************************************************************************
inline void CSamGenerator::store_mapping_result(
	uchar_t *id, uint32_t id_len, uint32_t flag, int32_t ref_seq_id, uint32_t ref_seq_pos, uint32_t mapping_quality, uchar_t* cigar, uint32_t *cigar_bin, uint32_t cigar_bin_len,
	int32_t mate_ref_seq_id, uint32_t mate_ref_seq_pos, int32_t template_len, genome_t dir, int32_t ref_length, uchar_t *sequence, uint32_t sequence_len,
	uchar_t *qualities, uint32_t edit_distance, uchar_t* mdz, double score, MatchingMethod method)
{
	if (store_BAM)
	{
		bam_part->SetRefId(ref_seq_id);
		if(ref_seq_id == -2)
			bam_part->SetPos(mate_ref_seq_pos);
		else
			bam_part->SetPos(ref_seq_pos);
		bam_part->SetRefLength(ref_length);
		bam_part->SetDirection(dir);

		bam_part->SetMapq(mapping_quality);
		bam_part->SetFlag(flag);
		bam_part->SetSeq(sequence, sequence_len);
		if (mate_ref_seq_id == -2)
			bam_part->SetNextRefId(ref_seq_id);
		else
			bam_part->SetNextRefId(mate_ref_seq_id);
		bam_part->SetNextPos(mate_ref_seq_pos);
		bam_part->SetTlen(template_len);
		bam_part->SetReadName(id+1);
		bam_part->SetCigar(cigar_bin, cigar_bin_len);
		bam_part->SetQual(qualities);

		bam_part->AddAuxInt(a_NM, edit_distance);
		if (mdz)
			bam_part->AddAuxString(a_MD, mdz);
		else
			bam_part->AddAuxString(a_MD, (uchar_t*) uc_star);

		if (params->read_group_id.length() > 0) {
			bam_part->AddAuxString(a_RG, (uchar_t*)params->read_group_id.c_str());
		}

#ifdef STORE_EXTRA_SAM_FIELDS
		bam_part->AddAuxFloat(a_AS, (float) score);

		string methodStr = method.toString();
		bam_part->AddAuxString(a_XM, methodStr);
#endif

		bam_part->CloseLine();
	}
	else
	{
		// ***************** Obligatory fields
		// Read id
		sam_part->AppendPart(id + 1, id_len - 1);

		// Flags
		sam_part->AppendPart(flag);

		// Ref. seq. name
		if (ref_seq_id == -1)
			sam_part->AppendPart("*");
		else if (ref_seq_id == -2)
			sam_part->AppendPart("=");
		else
			sam_part->AppendPart((*ref_seq_desc)[ref_seq_id].name);

		// Mapping position
		sam_part->AppendPart(ref_seq_pos);

		// Mapping quality
		sam_part->AppendPart(mapping_quality);

		// CIGAR
		if (cigar != nullptr) {
			sam_part->AppendPart(cigar);
		}
		else {
			sam_part->AppendPart("*");
		}

		// Mate ref. seq. name
		if (mate_ref_seq_id == -1)
			sam_part->AppendPart("*");
		else if (mate_ref_seq_id == -2 || mate_ref_seq_id == ref_seq_id)
			sam_part->AppendPart("=");
		else
			sam_part->AppendPart((*ref_seq_desc)[mate_ref_seq_id].name);

		// Mate position
		sam_part->AppendPart(mate_ref_seq_pos);

		// Template len
		sam_part->AppendPart(template_len);

		// Read sequence
		if (dir == genome_t::direct)
			sam_part->AppendPart(sequence, sequence_len);
		else
			sam_part->AppendPartRevComp(sequence, sequence_len);

		// Read qualities
		if (dir == genome_t::direct)
			sam_part->AppendPart(qualities, sequence_len);
		else
			sam_part->AppendPartRev(qualities, sequence_len);

		// ***************** Optional fields
		// Alignment score for mate (YS:i)

		// Edit distance (NM:i)
		sam_part->AppendPart("NM:i:", false);
		sam_part->AppendPart(edit_distance);

		// Match description (MD:Z)
		sam_part->AppendPart("MD:Z:", false);
		if (mdz)
			sam_part->AppendPart(mdz);
		else
			sam_part->AppendPart("*");

		if (params->read_group_id.length() > 0) {
			sam_part->AppendPart("RG:Z:", false);
			sam_part->AppendPart(params->read_group_id);
		}

#ifdef STORE_EXTRA_SAM_FIELDS
		sam_part->AppendPart("AS:i:", false);
		sam_part->AppendPart(score, 1);

		string methodStr = method.toString();
		sam_part->AppendPart("XM:Z:" + methodStr);

		sam_part->AppendPart("MO:Z:", false);
		sam_part->AppendPart(insertSizeModel->mean(), 13, false);
		sam_part->AppendPart("_", false);
		sam_part->AppendPart(insertSizeModel->dev(), 13);
#endif

		sam_part->CloseLine();
	}
}

// ************************************************************************************
void CSamGenerator::store_mapped_read(uchar_t *id, uchar_t *sequence, uchar_t *plus, uchar_t *quality,
	uint32_t id_len, uint32_t sequence_len, uint32_t plus_len, uint32_t quality_len,
	uchar_t *&data_ptr)
{
	uchar_t perfect_cigar[10] = { 0 };
	uchar_t perfect_mdz[10] = { 0 };
	uint32_t perfect_cigar_bin[1] = { 0 };
	uint32_t perfect_cigar_bin_len = 1;

	int end = CNumericConversions::Int2PChar(sequence_len, perfect_cigar);
	perfect_cigar[end] = 'M';
	perfect_cigar[end + 1] = 0;
	CNumericConversions::Int2PChar(sequence_len, perfect_mdz);
	perfect_mdz[end] = 0;
	perfect_cigar_bin[0] = sequence_len << 4;
	
	Evaluator<mapping_desc_t> hitEvaluator(*insertSizeModel);
	stored_mapped++;

	check_buffers_occupation();

	uint64_t no_mappings;
	LoadUInt(data_ptr+id_bytes, no_mappings, mapping_counter_size);

	data_ptr += (uint64_t) id_bytes + mapping_counter_size;

	uint64_t no_stored_mappings = no_mappings;
	//uint32_t rec_size = sizeof(ref_pos_t) + 1 + 1 + 1;

	string ref_seq_name;
	int32_t ref_seq_id;
	int32_t ref_seq_pos;
	string str;
	
	id_len = trim_id(id, id_len);

	mapping_desc[0].clear();

	bool heap_constructed = false;

	candidate_mapping_t candidate_mapping;

	for (uint64_t i = 0; i < no_stored_mappings; ++i)
	{
		ASSERT(data_ptr + candidate_mapping_t::max_size <= data_srt_end, "Error 2b");
		data_ptr += candidate_mapping.deserialize(data_ptr);

		if (ref_seq_desc->Translate(candidate_mapping.pos, ref_seq_name, ref_seq_pos, ref_seq_id, id))
		{
			uint32_t check_pos;
			ref_seq_desc->RevTranslate(check_pos, ref_seq_pos, ref_seq_id);

			//if (check_pos != candidate_mapping.pos) {}

			if (candidate_mapping.direction == genome_t::rev_comp)
				ref_seq_pos -= sequence_len - 1;
		}
		else
		{
			//cerr << "Error in position translation!\n";
			ref_seq_name = "??";
			ref_seq_pos = 0;
		}

		if (mapping_desc[0].size() == max_no_mappings && !heap_constructed)
		{
			make_heap(mapping_desc[0].begin(), mapping_desc[0].end(), mapping_desc_comparator);
			heap_constructed = true;
		}

		mapping_desc[0].push_back(mapping_desc_t(candidate_mapping, ref_seq_id, ref_seq_pos, sequence_len, nullptr,
			params->gap_ins_open, params->gap_ins_extend, params->gap_del_open, params->gap_del_extend, params->mismatch_score));
		mapping_desc[0].back().ref_length = sequence_len; // this will be corrected later
		mapping_desc[0].back().score = hitEvaluator(mapping_desc[0].back());

		if (heap_constructed)
		{
			push_heap(mapping_desc[0].begin(), mapping_desc[0].end(), mapping_desc_comparator);
			pop_heap(mapping_desc[0].begin(), mapping_desc[0].end(), mapping_desc_comparator);
			mapping_desc[0].pop_back();
		}

		ASSERT(data_ptr <= data_srt_end, "Error 3");
	}

	// fixme: this should be always true
	if (mapping_desc[0].size()) {

		std::vector<mapping_desc_t*> hits_se;

		// Sort mapping positions and remove redundant results 
		generate_unique_hits(mapping_desc[0], hits_se);

		MappingEvaluation se_eval(params->mapq_mult, params->mapq_div);
		std::stable_sort(hits_se.begin(), hits_se.end(), mapping_desc_t::compareByScoresAndPosDescending_ptr);
		se_eval.initialize(hits_se);

		mapping_desc_t& examined_mapping = *hits_se[0];

		if (examined_mapping.mapping.type == mapping_type_t::lev)
			fill_cigar_with_lev(examined_mapping, id, sequence, sequence_len, quality, false, true);
		else
			fill_cigar_with_indel_matching(examined_mapping, id, sequence, sequence_len, quality, true);

		softClipping->clipIllegalPositions(examined_mapping, insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);

		uchar_t* cigar = perfect_cigar;
		uchar_t* mdz = perfect_mdz;
		uint32_t *cigar_bin = perfect_cigar_bin;
		uint32_t cigar_bin_len = perfect_cigar_bin_len;

		if (examined_mapping.method != MatchingMethod::Enumeration::Perfect) {
			mp_mdz->Reserve(mdz);

			if (store_BAM)
			{
				mp_cigar_bin->Reserve(cigar_bin);
				cigar_bin_len = extract_cigarbin_and_mdz(examined_mapping.ext_cigar, cigar_bin, mdz);
			}
			else
			{
				mp_cigar->Reserve(cigar);
				extract_cigar_and_mdz(examined_mapping.ext_cigar, cigar, mdz);
			}
		}

		uint32_t flag = 0x0; // mapped single segment

		if (examined_mapping.mapping.direction == genome_t::rev_comp)
			flag += 0x10;							// read is reverse complemented
		
		store_mapping_result(
			id, id_len, flag, examined_mapping.ref_seq_id, examined_mapping.ref_seq_pos, examined_mapping.mapq, cigar, cigar_bin, cigar_bin_len,
			-1, 0, 0, examined_mapping.mapping.direction, examined_mapping.ref_length, sequence, sequence_len, quality,
			examined_mapping.err_edit_distance, mdz, examined_mapping.score, examined_mapping.method);

		if (examined_mapping.method != MatchingMethod::Enumeration::Perfect) {
			mp_mdz->Free(mdz);
			if(store_BAM)
				mp_cigar_bin->Free(cigar_bin);
			else
				mp_cigar->Free(cigar);
			mp_ext_cigar->Free(examined_mapping.ext_cigar);
		}
	}
	else {
		store_mapping_result(
			id, id_len, 0x04,
			-1, 0, 0, nullptr, nullptr, 0,
			-1, 0, 0,
			genome_t::direct, 1, sequence, sequence_len, quality,
			0, nullptr, 0, MatchingMethod::Enumeration::Unmatched);
	}	
}

// ************************************************************************************
bool  CSamGenerator::update_model(read_id_t read_id, uchar_t *id[2], uchar_t *sequence[2],
	uint32_t id_len[2], uint32_t sequence_len[2], uchar_t *&data_ptr)
{
	check_buffers_occupation();

	Evaluator<mapping_desc_t> hitEvaluator(*insertSizeModel);
	Evaluator<mapping_pair_t> pairEvaluator(*insertSizeModel);

	for (uint32_t r = 0; r < 2; ++r)
		id_len[r] = trim_id(id[r], id_len[r]);

	// Load mapping results for both reads
	for (uint32_t r = 0; r < 2; ++r)
	{
		mapping_desc[r].clear();
		uint64_t no_mappings;

		// Check whether there are any mapping results for current read id
		if (!find_id(read_id + r, data_ptr, data_srt_end))
			continue;

		LoadUInt(data_ptr + id_bytes, no_mappings, mapping_counter_size);
		ASSERT(data_ptr < data_srt_end, "Error 1");

		data_ptr += (uint64_t) id_bytes + mapping_counter_size;
		ASSERT(data_ptr < data_srt_end, "Error 2");

		uint64_t no_stored_mappings = no_mappings;

		int32_t ref_seq_id;
		string ref_seq_name;
		int32_t ref_seq_pos;

		bool heap_constructed = false;

		candidate_mapping_t candidate_mapping;

		for (uint64_t i = 0; i < no_stored_mappings; ++i)
		{
			ASSERT(data_ptr + candidate_mapping_t::max_size <= data_srt_end, "Error 2b");
			data_ptr += candidate_mapping.deserialize(data_ptr);

			if (ref_seq_desc->Translate(candidate_mapping.pos, ref_seq_name, ref_seq_pos, ref_seq_id, id[r]))
			{
				if (candidate_mapping.direction == genome_t::rev_comp)
					ref_seq_pos -= sequence_len[r] - 1;
			}
			else
			{
				//cerr << "Error in position translation!\n";
				ref_seq_name = "??";
				ref_seq_pos = 0;
			}

			if (mapping_desc[r].size() == max_no_mappings && !heap_constructed)
			{
				make_heap(mapping_desc[r].begin(), mapping_desc[r].end(), mapping_desc_comparator);
				heap_constructed = true;
			}

			mapping_desc[r].push_back(mapping_desc_t(candidate_mapping, ref_seq_id, ref_seq_pos, sequence_len[r], nullptr,
				params->gap_ins_open, params->gap_ins_extend, params->gap_del_open, params->gap_del_extend, params->mismatch_score));
			mapping_desc[r].back().ref_length = sequence_len[r]; // this will be corrected later
			mapping_desc[r].back().score = hitEvaluator(mapping_desc[r].back());

			if (heap_constructed)
			{
				push_heap(mapping_desc[r].begin(), mapping_desc[r].end(), mapping_desc_comparator);
				pop_heap(mapping_desc[r].begin(), mapping_desc[r].end(), mapping_desc_comparator);
				mapping_desc[r].pop_back();
			}

			ASSERT(data_ptr <= data_srt_end, "Error 3");
		}
	}

	std::vector<mapping_desc_t*> hits_se[2];

	// Sort mapping positions and remove redundant results 
	for (uint32_t i = 0; i < 2; ++i) {
		mapping_desc[i].reserve(2 * (mapping_desc[i].size() + mapping_desc[!i].size())); // for each hit we can find 2 new mate hits (myers, clipping) 
		generate_unique_hits(mapping_desc[i], hits_se[i]);

		/*
		auto it = find_if(hits_se[i].begin(), hits_se[i].end(), [](const mapping_desc_t* x)->bool {
			return x->err_edit_distance == 0;
		});
		*/

		hits_se[i].reserve(hits_se[i].size() + hits_se[!i].size());
	}

	// unique hits
	if (hits_se[0].size() == 1 && hits_se[1].size() == 1) {
		mapping_pair_t examined_mapping(hits_se[0][0], hits_se[1][0]);
		int tlen = examined_mapping.calculateInsertSize();
		if (tlen <= (int) max_mate_distance) {
			insertSizeModel->addSample(tlen);
		}
	}

	return insertSizeModel->getSamplesCount() < insertSizeModel->getMinSamplesCount();
}

// ************************************************************************************
void CSamGenerator::store_mapped_pair_reads(read_id_t read_id, uchar_t *id[2], uchar_t *sequence[2], uchar_t *plus[2], uchar_t *quality[2],
	uint32_t id_len[2], uint32_t sequence_len[2], uint32_t plus_len[2], uint32_t quality_len[2],
	uchar_t *&data_ptr)
{
	// **********************************
	// Initialisation
	CThreadWatch thr_watch;

	uchar_t perfect_cigar[2][10] = { { 0 },{ 0 } };
	uchar_t perfect_mdz[2][10] = { { 0 },{ 0 } };
	uint32_t perfect_cigar_bin[2][1] = { { 0 },{ 0 } };
	uint32_t perfect_cigar_bin_len[2] = { 1, 1 };

	for (int r = 0; r < 2; ++r) {
		int end = CNumericConversions::Int2PChar(sequence_len[r], perfect_cigar[r]);
		perfect_cigar[r][end] = 'M';
		perfect_cigar[r][end + 1] = 0;
		CNumericConversions::Int2PChar(sequence_len[r], perfect_mdz[r]);
		perfect_mdz[r][end] = 0;
		perfect_cigar_bin[r][0] = sequence_len[r] << 4;
	}
	
	check_buffers_occupation();

	Evaluator<mapping_desc_t> hitEvaluator(*insertSizeModel);
	Evaluator<mapping_pair_t> pairEvaluator(*insertSizeModel);

	for(uint32_t r = 0; r < 2; ++r)
		id_len[r] = trim_id(id[r], id_len[r]);

#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

	// **********************************
	// Load mapping results for both reads
	for(uint32_t r = 0; r < 2; ++r)
	{
		mapping_desc[r].clear();
		uint64_t no_mappings;

		// Check whether there are any mapping results for current read id
		if(!find_id(read_id+r, data_ptr, data_srt_end))
			continue;

		LoadUInt(data_ptr+id_bytes, no_mappings, mapping_counter_size);
		ASSERT(data_ptr < data_srt_end, "Error 1");
		
		data_ptr += (uint64_t) id_bytes + mapping_counter_size;
		ASSERT(data_ptr < data_srt_end, "Error 2");
	
		int32_t ref_seq_id;
		string ref_seq_name;
		int32_t ref_seq_pos;

		bool heap_constructed = false;
		
		candidate_mapping_t candidate_mapping;

		for(uint64_t i = 0; i < no_mappings; ++i)
		{
			ASSERT(data_ptr + candidate_mapping_t::max_size <= data_srt_end, "Error 2b");
			data_ptr += candidate_mapping.deserialize(data_ptr);

			if(ref_seq_desc->Translate(candidate_mapping.pos, ref_seq_name, ref_seq_pos, ref_seq_id, id[r]))
			{
				if(candidate_mapping.direction == genome_t::rev_comp)
					ref_seq_pos -= sequence_len[r] - 1;
			}
			else
			{
				//cerr << "Error in position translation!\n";
				ref_seq_name = "??";
				ref_seq_pos = 0;
			}

			// if lots of mapping (max_no_mappings limit reached) - use heap to select best
			if (mapping_desc[r].size() == max_no_mappings && !heap_constructed)
			{
				make_heap(mapping_desc[r].begin(), mapping_desc[r].end(), mapping_desc_comparator);
				heap_constructed = true;
			}

			mapping_desc[r].push_back(mapping_desc_t(candidate_mapping, ref_seq_id, ref_seq_pos, sequence_len[r], nullptr, 
				params->gap_ins_open, params->gap_ins_extend, params->gap_del_open, params->gap_del_extend, params->mismatch_score));
			mapping_desc[r].back().ref_length = sequence_len[r]; // this will be corrected later
			mapping_desc[r].back().score = hitEvaluator(mapping_desc[r].back());

			if (heap_constructed)
			{
				push_heap(mapping_desc[r].begin(), mapping_desc[r].end(), mapping_desc_comparator);
				pop_heap(mapping_desc[r].begin(), mapping_desc[r].end(), mapping_desc_comparator);
				mapping_desc[r].pop_back();
			}

			ASSERT(data_ptr <= data_srt_end, "Error 3");
		}
	}

	// Remove mappings with too much errors for short reads
	for(int r = 0; r < 2; ++r)
		if (sequence_len[r] < params->min_read_len)
		{
			if (!params->enable_short_reads)// || (sequence_len[!r] != sequence_len[r]))
				mapping_desc[r].clear();
			else
			{
				uint32_t err_thr = (uint32_t)(params->max_no_errors * params->sensitivity_factor * sequence_len[r] / params->min_read_len);

				auto new_end = remove_if(mapping_desc[r].begin(), mapping_desc[r].end(), [err_thr](mapping_desc_t md) {
					return md.err_edit_distance > err_thr;
					});

				mapping_desc[r].erase(new_end, mapping_desc[r].end());
			}
		}

#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[0] += thr_watch.GetElapsedTime();
#endif

	// **********************************
	// Generation of unique hits
#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

	uint32_t exact_counter = 0;
	std::vector<mapping_desc_t*> hits_se[2];

	// Sort mapping positions and remove redundant results 
	for(uint32_t i = 0; i < 2; ++i) {
		mapping_desc[i].reserve(16 * (mapping_desc[i].size() + mapping_desc[!i].size())); // for each hit we can find 2 new mate hits (myers, clipping) 
		generate_unique_hits(mapping_desc[i], hits_se[i]);
		
		auto it = find_if(hits_se[i].begin(), hits_se[i].end(), [](const mapping_desc_t* x)->bool {
			return x->err_edit_distance == 0;
		});

		if (it != hits_se[i].end())
			++exact_counter;

		hits_se[i].reserve(hits_se[i].size() + hits_se[!i].size());
	}

#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[1] += thr_watch.GetElapsedTime();
#endif
 
	// **********************************
#	// Matching mappings obtained in the mapping stages
	// find pairings
	std::pair<int, int> entireRange(0, std::numeric_limits<int>::max());
	std::pair<int, int> lowErrorRange(0, params->max_no_errors);
	std::pair<int, int> highErrorRange(params->max_no_errors + 1, std::numeric_limits<int>::max());

	// Try find by matching pairs from mapping stage
	std::vector<mapping_pair_t> close_mapping_pairs;
	std::vector<mapping_pair_t> variant_mapping_pairs;
	std::vector<mapping_pair_t> myers_mapping_pairs;
	std::vector<mapping_pair_t> clipping_mapping_pairs;

#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

	bool close_found = find_pairs<MatchingMethod::Enumeration::FromMapping>(mapping_desc, hits_se, id, sequence, sequence_len, quality, close_mapping_pairs, entireRange);
		
	bool ignore_pair = false;
	if (!close_found && mapping_desc[0].size() > params->max_no_mappings / 2 && mapping_desc[1].size() > params->max_no_mappings / 2)
		ignore_pair = true;

#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[2] += thr_watch.GetElapsedTime();
#endif

	// **********************************
	// Looking for matches using external variant database
#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

#ifdef ENABLE_VCF_VARIANTS
	if (!close_found && !ignore_pair && params->enable_var_indel_long)
	{
		std::vector<mapping_pair_t> variant_mapping_pairs_ins;
		std::vector<mapping_pair_t> variant_mapping_pairs_del;

		close_found |= find_pairs<MatchingMethod::Enumeration::VariantInsertLong>(mapping_desc, hits_se, id, sequence, sequence_len, quality, variant_mapping_pairs_ins, lowErrorRange);
		close_found |= find_pairs<MatchingMethod::Enumeration::VariantDeleteLong>(mapping_desc, hits_se, id, sequence, sequence_len, quality, variant_mapping_pairs_del, lowErrorRange);

		variant_mapping_pairs.resize(variant_mapping_pairs_ins.size() + variant_mapping_pairs_del.size());
		std::merge(variant_mapping_pairs_ins.begin(), variant_mapping_pairs_ins.end(),
			variant_mapping_pairs_del.begin(), variant_mapping_pairs_del.end(),
			variant_mapping_pairs.begin(), mapping_pair_t::compareByScoresDescending);

		if (close_found)
			close_mapping_pairs = std::move(variant_mapping_pairs);
	}
#endif

#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[3] += thr_watch.GetElapsedTime();
#endif

	// **********************************
	// Looking for matches using paired Myers and/or clipping - low error range
	// Try rescue pairs by LevMyers and clipping search in mapping regions
	if (!close_found && !ignore_pair) {
#ifdef MEASURE_PP_TIMES
		thr_watch.StartTimer();
#endif

		close_found |= find_pairs<MatchingMethod::Enumeration::LevMyers>(mapping_desc, hits_se, id, sequence, sequence_len, quality, myers_mapping_pairs, lowErrorRange);

		bool tryClipping = true;
		if (myers_mapping_pairs.size()) {
			const auto &best = myers_mapping_pairs.front();
			if ((best.first->method == MatchingMethod::Enumeration::LevMyers && best.first->err_edit_distance <= 1) ||
				(best.second->method == MatchingMethod::Enumeration::LevMyers && best.second->err_edit_distance <= 1)) {
				tryClipping = false;
			}
		}

#ifdef MEASURE_PP_TIMES
		thr_watch.StopTimer();
		part_times[4] += thr_watch.GetElapsedTime();
#endif

		if (tryClipping) {
#ifdef MEASURE_PP_TIMES
			thr_watch.StartTimer();
#endif

			close_found |= find_pairs<MatchingMethod::Enumeration::Clipping>(mapping_desc, hits_se, id, sequence, sequence_len, quality, clipping_mapping_pairs, lowErrorRange);
			
			close_mapping_pairs.resize(myers_mapping_pairs.size() + clipping_mapping_pairs.size());
			std::merge(myers_mapping_pairs.begin(), myers_mapping_pairs.end(),
				clipping_mapping_pairs.begin(), clipping_mapping_pairs.end(),
				close_mapping_pairs.begin(), mapping_pair_t::compareByScoresDescending);

#ifdef MEASURE_PP_TIMES
			thr_watch.StopTimer();
			part_times[5] += thr_watch.GetElapsedTime();
#endif

		}
		else {
			close_mapping_pairs = std::move(myers_mapping_pairs);
		}
	}

	// **********************************
	// Looking for matches using paired Myers and/or clipping - high error range
	if (!close_found && !ignore_pair) {
#ifdef MEASURE_PP_TIMES
		thr_watch.StartTimer();
#endif
		close_found |= find_pairs<MatchingMethod::Enumeration::LevMyers>(mapping_desc, hits_se, id, sequence, sequence_len, quality, myers_mapping_pairs, highErrorRange);
#ifdef MEASURE_PP_TIMES
		thr_watch.StopTimer();
		part_times[6] += thr_watch.GetElapsedTime();
#endif
		
#ifdef MEASURE_PP_TIMES
		thr_watch.StartTimer();
#endif
		close_found |= find_pairs<MatchingMethod::Enumeration::Clipping>(mapping_desc, hits_se, id, sequence, sequence_len, quality, clipping_mapping_pairs, highErrorRange);
#ifdef MEASURE_PP_TIMES
		thr_watch.StopTimer();
		part_times[7] += thr_watch.GetElapsedTime();
#endif

		close_mapping_pairs.resize(myers_mapping_pairs.size() + clipping_mapping_pairs.size());

		std::merge(myers_mapping_pairs.begin(), myers_mapping_pairs.end(),
			clipping_mapping_pairs.begin(), clipping_mapping_pairs.end(),
			close_mapping_pairs.begin(), mapping_pair_t::compareByScoresDescending);
	}

	// **********************************
	// Checking global best SE mappings as a pair
	// find best globally and check if it is close
	if (!close_found && hits_se[0].size() && hits_se[1].size()) {
		mapping_desc_t* bestHits[2];
		bestHits[0] = *std::min_element(hits_se[0].begin(), hits_se[0].end(), mapping_desc_t::compareByScoresAndPosDescending_ptr);
		bestHits[1] = *std::min_element(hits_se[1].begin(), hits_se[1].end(), mapping_desc_t::compareByScoresAndPosDescending_ptr);

		mapping_pair_t mapping(bestHits[0], bestHits[1]);
		mapping.score = pairEvaluator(mapping);

		if (mapping.calculateInsertSize() <= insertSizeModel->mean() + insertSizeModel->getMyersDev()) {
			close_mapping_pairs.push_back(mapping);
			close_found = true;
		}
	}

	// **********************************
	// Looking for distant matches
	// Check whether distant mapping is better than hits already found 
	std::vector<mapping_pair_t> distant_mapping_pairs;
	bool distant_found = false;

#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

	for(int r = 0; r < 2; ++r)
		std::stable_sort(hits_se[r].begin(), hits_se[r].end(), mapping_desc_t::compareByPositions_ptr);

	distant_found = find_pairs<MatchingMethod::Enumeration::FromMappingDistant>(mapping_desc, hits_se, id, sequence, sequence_len, quality, distant_mapping_pairs, entireRange);
#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[8] += thr_watch.GetElapsedTime();
#endif

	// **********************************
	// Calculation of SE qualities
#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

	// calculate single-end qualities
	std::vector<MappingEvaluation> se_eval;
	se_eval.emplace_back(params->mapq_mult, params->mapq_div);
	se_eval.emplace_back(params->mapq_mult, params->mapq_div);
	for (int i = 0; i < 2; ++i) {
		if (hits_se[i].size()) {
			std::stable_sort(hits_se[i].begin(), hits_se[i].end(), mapping_desc_t::compareByScoresAndPosDescending_ptr);
			se_eval[i].initialize(hits_se[i]);	
		}
	}
#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[9] += thr_watch.GetElapsedTime();
#endif

	// If some paired-end mappings were found
	MappingEvaluation pe_eval(params->mapq_mult, params->mapq_div);
	std::vector<mapping_pair_t> mapping_pairs;
	int pair_id = -1;
	if (close_found || distant_found) {

		// **********************************
		// Merging found paired mapping
#ifdef MEASURE_PP_TIMES
		thr_watch.StartTimer();
#endif

		// Merge close and distant mappings
		if (!close_found) {
			mapping_pairs = std::move(distant_mapping_pairs);
		} else if (!distant_found) {
			mapping_pairs = std::move(close_mapping_pairs);
		} else {
			mapping_pairs.resize(close_mapping_pairs.size() + distant_mapping_pairs.size());
			std::merge(close_mapping_pairs.begin(), close_mapping_pairs.end(),
				distant_mapping_pairs.begin(), distant_mapping_pairs.end(),
				mapping_pairs.begin(), mapping_pair_t::compareByScoresDescending);
		}
	
		sort(mapping_pairs.begin(), mapping_pairs.end(), mapping_pair_t::compareByScoresAndPosDescending);

		//find_duplicate_pairs(mapping_pairs);

		// calculate paired-end qualities
		pe_eval.initialize(*insertSizeModel, params->score_discretization_threshold, mapping_pairs);

		// get range of equally good mappings
		std::uniform_int_distribution<int> distribution(0, pe_eval.best[0].count - 1);
		pair_id = distribution(engine);
		
		// perform clipping of illegal positions	
		auto&mp = mapping_pairs[pair_id];
		for (int r = 0; r < 2; ++r) {
			// remove read contribution from pair score
			mp.score -= mp[r].score;
			
			if (mp[r].method == MatchingMethod::Enumeration::Unmatched
				|| !mp[r].ext_cigar || mp[r].ext_cigar[0] == 0)
			{
				int isize = mp.calculateInsertSize();
				bool isDistant = 
					isize > insertSizeModel->mean() + insertSizeModel->getMyersDev() ||
					isize < insertSizeModel->mean() - insertSizeModel->getMyersDev();

				fill_cigar(mp[r], id[r], sequence[r], sequence_len[r], quality[r], isDistant, true);
			}

			softClipping->clipIllegalPositions(mp[r], insertSizeModel->getScoring(), ref_seq_desc->GetDescription(), *mp_ext_cigar);
			if (params->enable_boundary_clipping) {
				softClipping->clipBoundaryPositions(mp[r], insertSizeModel->getScoring());
			}

			// restore read contribution to pair score
			mp.score += mp[r].score;
		}
		
		if (params->enable_paired_clipping) {
			softClipping->clipOverlappingPairs(mp, insertSizeModel->getScoring(), *mp_ext_cigar);
		}

#ifdef MEASURE_PP_TIMES
		thr_watch.StopTimer();
		part_times[10] += thr_watch.GetElapsedTime();
#endif

	}
	// find single mappings if everything other fails	
	else  {
		// ********************************** 
		// Finding best SE mappings
#ifdef MEASURE_PP_TIMES
		thr_watch.StartTimer();
#endif
		find_single(hits_se, id, sequence, sequence_len, quality);
#ifdef MEASURE_PP_TIMES
		thr_watch.StopTimer();
		part_times[11] += thr_watch.GetElapsedTime();
#endif

	}

	bool filterPassed = params->filter.length() == 0; // if filter not specified - assume it is passed
	//bool mappings_removed = false;

#ifdef MEASURE_PP_TIMES
	thr_watch.StartTimer();
#endif

	// Store mapping results for pairs if there is at least one pair
	if (mapping_pairs.size())
	{
		// **********************************
		// Storing paired mapping
		stored_mapped += 2;
		bool first_pair = true;

		auto &p = mapping_pairs[pair_id]; // take into account only first mapping
		{
			++stat_mapped_pe_method[(int)p.first->method];
			++stat_mapped_pe_method[(int)p.second->method];

			int32_t template_len = p.calculateInsertSize();

			// use unique mappings to update insert size model
			if (
				(p.first->method == MatchingMethod::Enumeration::FromMapping || p.first->method == MatchingMethod::Enumeration::Perfect) &&
				(p.second->method == MatchingMethod::Enumeration::FromMapping || p.second->method == MatchingMethod::Enumeration::Perfect) &&
				template_len != mapping_pair_t::MAX_INSERT_SIZE &&
				pe_eval.best[0].count == 1 && pe_eval.best[1].count == 0) {
			
				if (template_len <= (int)max_mate_distance) {
					//					cout << "tlen: " + to_string(template_len) + "\n";
					insertSizeModel->addSample(template_len);
				}
			}

			int rnext[2] = { -2, -2 };

			// set flags
			uint32_t flag[2] = { first_pair ? 0x3u : 0x103u, first_pair ? 0x3u : 0x103u };	// both paired reads mapped; 0x100 - not the first mapping (more tha one mapping found)

			if (template_len == mapping_pair_t::MAX_INSERT_SIZE) {
				template_len = 0;
				rnext[0] = p.second->ref_seq_id;
				rnext[1] = p.first->ref_seq_id;
				flag[0] = flag[1] = first_pair ? 0x1 : 0x101; // if mapped to different chromosomes 
			}

			for (int r = 0; r < 2; ++r) {
				auto& examined_mapping = p[r];

				if (examined_mapping.mapping.direction == genome_t::rev_comp) {
					flag[r] += 0x10;
					flag[!r] += 0x20;
				}
				flag[r] += (r == 0) ? 0x40 : 0x80;		// first/second read of a pair

				if (examined_mapping.err_edit_distance == 0) { --exact_counter; }

				// if chromosome name fits filter
				if ((*ref_seq_desc)[examined_mapping.ref_seq_id].name == params->filter) {
					filterPassed |= true;
				}
			}

			for (int r = 0; r < 2; ++r) {
				auto& examined_mapping = p[r];
				auto& mate_mapping = p[!r];
				int insert_size = (examined_mapping.mapping.pos < mate_mapping.mapping.pos) ? template_len : -template_len;

				uchar_t* cigar = perfect_cigar[r];
				uchar_t* mdz = perfect_mdz[r];
				uint32_t* cigar_bin = perfect_cigar_bin[r];
				uint32_t cigar_bin_len = perfect_cigar_bin_len[r];

				if (examined_mapping.method != MatchingMethod::Enumeration::Perfect && examined_mapping.method != MatchingMethod::Enumeration::PerfectDistant) {
					mp_mdz->Reserve(mdz);
					if (store_BAM)
					{
						mp_cigar_bin->Reserve(cigar_bin);
						cigar_bin_len = extract_cigarbin_and_mdz(examined_mapping.ext_cigar, cigar_bin, mdz);
					}
					else
					{
						mp_cigar->Reserve(cigar);
						extract_cigar_and_mdz(examined_mapping.ext_cigar, cigar, mdz);
					}

#ifdef REFINE_SCORE_SNP
					refine_score_snp(examined_mapping, id[r], id_len[r], sequence[r], sequence_len[r], examined_mapping.mapping.direction);
#endif
#ifdef VALIDATE_CIGAR
					validate_cigar(examined_mapping, id[r], id_len[r], sequence[r], sequence_len[r], examined_mapping.mapping.direction);
#endif
				}

				if (filterPassed) {
					store_mapping_result(
						id[r], id_len[r], flag[r], examined_mapping.ref_seq_id, examined_mapping.ref_seq_pos, examined_mapping.mapq, cigar, cigar_bin, cigar_bin_len,
						rnext[r], mate_mapping.ref_seq_pos, insert_size, examined_mapping.mapping.direction, examined_mapping.ref_length, sequence[r], sequence_len[r], quality[r],
						examined_mapping.err_edit_distance, mdz, examined_mapping.score, examined_mapping.method);
				}

				if (examined_mapping.method != MatchingMethod::Enumeration::Perfect && examined_mapping.method != MatchingMethod::Enumeration::PerfectDistant) {
					mp_mdz->Free(mdz);
					if(store_BAM)
						mp_cigar_bin->Free(cigar_bin);
					else
						mp_cigar->Free(cigar);
				}
			}

#if 0
			if (mapping_pairs.size() > 100)
			{
				cout << string((char*)(id[0])) + "  :  " +
					"MP_size: " + to_string(mapping_pairs.size()) +
					"   Hits1:" + to_string(mapping_desc[0].size()) +
					"   Hits2:" + to_string(mapping_desc[1].size()) +
					"   MAPQ1: " + to_string(p[0].mapq) +
					"   MAPQ2: " + to_string(p[1].mapq) + 
					"   Meth1: " + MatchingMethod(mapping_pairs.front()[0].method).toString() +
					"   Meth2: " + MatchingMethod(mapping_pairs.front()[1].method).toString() +
					"   E-dist1: " + to_string(mapping_desc[0].front().err_edit_distance) +
					"   E-dist2: " + to_string(mapping_desc[1].front().err_edit_distance) +
					"   E-dist p1: " + to_string(mapping_pairs.front()[0].err_edit_distance) +
					"   E-dist p2: " + to_string(mapping_pairs.front()[1].err_edit_distance) +
					"   Score0: " + to_string(mapping_pairs[0].score) +
					"   Score1: " + to_string(mapping_pairs[1].score) +
					"   PE-eval0: " + to_string(pe_eval.best[0].score) + " (" + to_string(pe_eval.best[0].count) + ")" + (mapping_pairs.front().is_distant ? "D" : "-") +
					"   PE-eval1: " + to_string(pe_eval.best[1].score) + " (" + to_string(pe_eval.best[1].count) + ")" +
					"\n";
				for (int i = 0; i < 2; ++i)
				{
					sort(hits_se[i].begin(), hits_se[i].end(), mapping_desc_t::compareByPositions_ptr);
					cout << "*** ";
					cout << to_string(hits_se[i].front()->raw_pos) + " ";

					for (int j = 1; j < hits_se[i].size(); ++j)
						cout << to_string(hits_se[i][j]->raw_pos - hits_se[i][j - 1]->raw_pos) + " (" +
						to_string(hits_se[i][j]->ref_seq_id) + " " + to_string(hits_se[i][j]->ref_seq_pos) + ")  ";
					cout << endl;
					sort(hits_se[i].begin(), hits_se[i].end(), mapping_desc_t::compareByScoresAndPosDescending_ptr);
				}
			}
#endif

			first_pair = false;
		}

		stat_mapped_pe_histo[stat_mapped_pe_histo.size() - 1] += exact_counter;

		if (pe_eval.best[0].count == 1)
			stat_mapped_pe_pair_unq++;
		else
			stat_mapped_pe_pair_more++;
	}
	else
	{
		// **********************************
		// Storing SE mappings
		// Store reads as mapped but unpaired
		stored_unmapped += 2;

		// check filter
		for (uint32_t r = 0; r < 2; ++r) {
			auto& mapping = hits_se[r];
			if (mapping.size() > 0 && (*ref_seq_desc)[mapping[0]->ref_seq_id].name == params->filter) {
				filterPassed |= true;
			}
		}

		// Log unmapped pairs
#ifdef LOG_UNMAPPED
		if (hits_se[0].empty() && hits_se[1].empty())
		{
			for (int i = 0; i < 2; ++i)
			{
				fwrite(id[i], 1, id_len[i], f_out[i]);
				putc('\n', f_out[i]);
				fwrite(sequence[i], 1, sequence_len[i], f_out[i]);
				putc('\n', f_out[i]);
				putc('+', f_out[i]);
				putc('\n', f_out[i]);
				fwrite(quality[i], 1, sequence_len[i], f_out[i]);
				putc('\n', f_out[i]);
			}
		}
#endif

		// write only when filter passed
		for(uint32_t r = 0; r < 2; ++r)
		{
			auto& current_mapping = hits_se[r];
			auto& mate_mapping = hits_se[!r];
			
			bool first = true;
						
			if (current_mapping.size() == 0)		// Read has no alignments
			{
				uint32_t flag = 0x01;		// the read is one of a pair
				flag += 0x04;				// the read is unmapped
				flag += r == 0 ? 0x40 : 0x80;				// first or last read of a pair

				if (mate_mapping.size() == 0) {
					// if mate is also unmapped
					flag += 0x08;
					
					if (filterPassed) {
						store_mapping_result(
							id[r], id_len[r], flag,
							-1, 0, 0, nullptr, nullptr, 0,
							-1, 0, 0,
							genome_t::direct, 1, sequence[r], sequence_len[r], quality[r],
							0, nullptr, 0, MatchingMethod::Enumeration::Unmatched);
					}
				}
				else  {
					// if mate is mapped - use RNAME and POS from mate (SAM recommended practices)
					if (filterPassed) {
						store_mapping_result(
							id[r], id_len[r], flag,
							mate_mapping[0]->ref_seq_id, mate_mapping[0]->ref_seq_pos, 0, nullptr, nullptr, 0,
							mate_mapping[0]->ref_seq_id, mate_mapping[0]->ref_seq_pos, 0,
							genome_t::direct, mate_mapping[0]->ref_length, sequence[r], sequence_len[r], quality[r],
							0, nullptr, 0, MatchingMethod::Enumeration::Unmatched);
					}
				}

				continue;
			}			

			for(auto &hit : current_mapping)
			{
				auto& examined_mapping = *hit;
				int32_t template_len = 0;
			
				uint32_t flag = first ? 0x1 : 0x101;		// one read of a pair; 0x100 - there is more than one mapping and the stored is not the first one
				flag += (examined_mapping.mapping.direction == genome_t::rev_comp) ? 0x10 : 0; // add 0x10 if read is reverse complemented					
				flag += (mate_mapping.size() == 0) ? 0x08 : 0;			// add 0x08 if the other read of the pair is unmapped
				flag += (r == 0) ? 0x40 : 0x80;							// first or last read of a pair
			
				uchar_t* cigar = perfect_cigar[r];
				uchar_t* mdz = perfect_mdz[r];
				uint32_t* cigar_bin = perfect_cigar_bin[r];
				uint32_t cigar_bin_len = perfect_cigar_bin_len[r];

				if (examined_mapping.method != MatchingMethod::Enumeration::Perfect && examined_mapping.method != MatchingMethod::Enumeration::PerfectDistant) {
					mp_mdz->Reserve(mdz);
					if (store_BAM)
					{
						mp_cigar_bin->Reserve(cigar_bin);
						cigar_bin_len = extract_cigarbin_and_mdz(examined_mapping.ext_cigar, cigar_bin, mdz);
					}
					else
					{
						mp_cigar->Reserve(cigar);
						extract_cigar_and_mdz(examined_mapping.ext_cigar, cigar, mdz);
					}
				}

				if (filterPassed) {
					store_mapping_result(
						id[r], id_len[r], flag, examined_mapping.ref_seq_id, examined_mapping.ref_seq_pos, examined_mapping.mapq, cigar, cigar_bin, cigar_bin_len,
						-1, 0, template_len, examined_mapping.mapping.direction, examined_mapping.ref_length, sequence[r], sequence_len[r], quality[r],
						examined_mapping.err_edit_distance, mdz, examined_mapping.score, examined_mapping.method);
				}

				if (examined_mapping.method != MatchingMethod::Enumeration::Perfect && examined_mapping.method != MatchingMethod::Enumeration::PerfectDistant) {
					mp_mdz->Free(mdz);
					if(store_BAM)
						mp_cigar_bin->Free(cigar_bin);
					else
						mp_cigar->Free(cigar);
				}

				first = false;
			}
		}

		if(hits_se[0].empty() && hits_se[1].empty())
			stat_mapped_pe_independent_none++;
		else if(hits_se[0].empty() || hits_se[1].empty())
			stat_mapped_pe_independent_single++;
		else
			stat_mapped_pe_independent_both++;
	}

#ifdef MEASURE_PP_TIMES
	thr_watch.StopTimer();
	part_times[12] += thr_watch.GetElapsedTime();
#endif
	// **********************************
	// Releasing cigars
	// clear all EXT-CIGARS
	for (int i = 0; i < 2; ++i)
		for (auto it = mapping_desc[i].begin(); it != mapping_desc[i].end(); ++it) {
			auto& hit = *it;
			if (hit.ext_cigar != nullptr) {
				mp_ext_cigar->FreeConditional(hit.ext_cigar);
				hit.ext_cigar = nullptr; // to prevent multiple deallocations
			}
		}
}

// ************************************************************************************
void CSamGenerator::store_unmapped_read(uchar_t *id, uchar_t *sequence, uchar_t *plus, uchar_t *quality,
	uint32_t id_len, uint32_t sequence_len, uint32_t plus_len, uint32_t quality_len,
	uchar_t *&data_ptr)
{
	stored_unmapped++;

	check_buffers_occupation();

	if (store_BAM)
	{
		bam_part->SetRefId(-1);
		bam_part->SetRefId(-1);
		bam_part->SetDirection(genome_t::direct);
		bam_part->SetMapq(0);
		bam_part->SetCigar(nullptr, 0);
		bam_part->SetFlag(4);		// flag
		bam_part->SetNextRefId(-1);
		bam_part->SetNextPos(-1);
		bam_part->SetTlen(0);

		bam_part->SetReadName(id+1);
		bam_part->SetSeq(sequence, sequence_len);
		bam_part->SetQual(quality);

		bam_part->CloseLine();
	}
	else
	{
		// ***************** Obligatory fields
		// Read id
		id_len = trim_id(id, id_len);
		sam_part->AppendPart(id + 1, id_len - 1);

		// Flags - unmapped read = 4
		sam_part->AppendPart(4);

		// Ref. seq. name
		sam_part->AppendPart("*");

		// Mapping position
		sam_part->AppendPart(0);

		// Mapping quality
		sam_part->AppendPart(0);

		// CIGAR
		sam_part->AppendPart("*");

		// Mate ref. seq. name
		sam_part->AppendPart("*");

		// Mate position
		sam_part->AppendPart(0);

		// Template len
		sam_part->AppendPart(0);

		// Read sequence
		sam_part->AppendPart(sequence, sequence_len);

		// Read qualities
		sam_part->AppendPart(quality, quality_len);

		sam_part->CloseLine();
	}
}
// ************************************************************************************
// Convert sequence to reverse complement
void CSamGenerator::convert_to_rev_comp(uchar_t *dest, uchar_t *src, uint32_t len)
{
	uint32_t pos = 0;

	for(int32_t i = len-1; i >= 0; --i)
	{
		switch(src[i])
		{
		case 'A':
			dest[pos++] = sym_code_T;
			break;
		case 'C':
			dest[pos++] = sym_code_G;
			break;
		case 'G':
			dest[pos++] = sym_code_C;
			break;
		case 'T':
			dest[pos++] = sym_code_A;
			break;
		default:
			dest[pos++] = sym_code_N_read;
			break;
		}
	}
}

// ************************************************************************************
void CSamGenerator::copy_direct(uchar_t *dest, uchar_t *src, uint32_t len)
{
	for(int32_t i = 0; i < (int32_t) len; )
	{
		switch(src[i])
		{
		case 'A':
			dest[i++] = sym_code_A;
			break;
		case 'C':
			dest[i++] = sym_code_C;
			break;
		case 'G':
			dest[i++] = sym_code_G;
			break;
		case 'T':
			dest[i++] = sym_code_T;
			break;
		default:
			dest[i++] = sym_code_N_read;
			break;
		}
	}
}

// ************************************************************************************
void CSamGenerator::ref_copy_direct(uchar_t *dest, uchar_t *src, ref_pos_t pos, uint32_t len)
{
	ref_pos_t packed_pos = pos / 2;
	uchar_t *packed_src = src + packed_pos;
	uchar_t* dest0 = dest;
	uint32_t len0 = len;

	if(pos & 1)
	{
		*dest++ = *packed_src++ & 0x0f;
		--len;
	}

	for(uint32_t i = 0; i < len / 2; ++i)
	{
		*dest++ = *packed_src >> 4;
		*dest++ = *packed_src++ & 0x0f;
	}

	if(len & 1)
		*dest = *packed_src >> 4;

	for (uint32_t i = 0; i < len0; ++i)
		if (dest0[i] > sym_code_N_ref)
			dest0[i] = sym_code_N_ref;
}

// ************************************************************************************
void CSamGenerator::check_buffers_occupation()
{
	if (store_BAM)
	{
		if (bam_part->IsFilled())
		{
			push_mapped_part();
			bam_part->Reserve();
		}
	}
	else
	{
		if (sam_part->IsFilled())
		{
			push_mapped_part();
			sam_part->Reserve();
		}
	}
}

#ifdef VALIDATE_CIGAR
// ************************************************************************************
bool CSamGenerator::validate_cigar(mapping_desc_t& md, uchar_t* id, uint32_t id_len, uchar_t* sequence, uint32_t sequence_len, genome_t direction)
{
	auto ref = reference->GetData();

	vector<uchar_t> genome;
	vector<uchar_t> read;

	uint32_t begin;
	uint32_t size = 2 * sequence_len + 2 * params->max_approx_indel_len + 2;
	genome.resize(size + 2ll);

	ref_seq_desc->RevTranslate(begin, md.ref_seq_pos, md.ref_seq_id);

	// Genome prefetch
	uint32_t i_genome = 0;
	uchar_t* genome_ptr = ref + (begin >> 1);
	uint32_t text_size_div2 = size / 2;

	if (begin & 1)
		genome[i_genome++] = "ACGTNNNNNNNNNNNN"[*genome_ptr++ & 0x0f];
	for (uint32_t i = 0; i < text_size_div2; ++i)
	{
		genome[i_genome++] = "ACGTNNNNNNNNNNNN"[*genome_ptr >> 4];
		genome[i_genome++] = "ACGTNNNNNNNNNNNN"[*genome_ptr++ & 0x0f];
	}

	read.resize(sequence_len);
	if (direction == genome_t::direct)
		for (int i = 0; i < sequence_len; ++i)
			read[i] = sequence[i];
	else
		for (int i = 0; i < sequence_len; ++i)
		{
			switch (sequence[i])
			{
			case 'A':
				read[sequence_len - 1 - i] = 'T'; break;
			case 'C':
				read[sequence_len - 1 - i] = 'G'; break;
			case 'G':
				read[sequence_len - 1 - i] = 'C'; break;
			case 'T':
				read[sequence_len - 1 - i] = 'A'; break;
			default:
				read[sequence_len - 1 - i] = 'N';
			}
		}

	bool is_correct = true;
	i_genome = 0;
	uint32_t i_read = 0;
	uchar_t* ec = md.ext_cigar;

	while(*ec)
	{
		if (*ec == '.')
		{
			if (genome[i_genome] != read[i_read])
			{
				is_correct = false;
				break;
			}
			++i_genome;
			++i_read;
			++ec;
		}
		else if (*ec == 'A' || *ec == 'C' || *ec == 'G' || *ec == 'T' || *ec == 'N')
		{
			if (genome[i_genome] != *ec || (genome[i_genome] == read[i_read] && *ec != 'N'))
			{
				is_correct = false;
				break;
			}
			++i_genome;
			++i_read;
			++ec;
		}
		else if (*ec == '$')
		{
			++i_read;
			++ec;
		}
		else if (*ec == '#')
		{
			++ec;
			if (!*ec)
			{
				is_correct = false;
				break;
			}

			if (read[i_read] != *ec)
			{
				is_correct = false;
				break;
			}
			++i_read;
			++ec;
		}
		else if (*ec == '^')
		{
			++ec;
			if (!*ec)
			{
				is_correct = false;
				break;
			}

			if (genome[i_genome] != *ec)
			{
				is_correct = false;
				break;
			}
			++i_genome;
			++ec;
		}
	}

	if (!is_correct || i_read != sequence_len)
	{
		string err;
		
		if(i_read != sequence_len)
			err = "Error in ExtCigar (wrong length): " + string(id, id + id_len) + " :  chr: " + to_string(md.ref_seq_id) + " pos: " + to_string(md.ref_seq_pos) + "  ec: ";
		else
			err = "Error in ExtCigar: " + string(id, id + id_len) + " :  chr: " + to_string(md.ref_seq_id) + " pos: " + to_string(md.ref_seq_pos) + "  ec: ";
		err.push_back(*ec);
		if (md.mapping.direction == genome_t::direct)
			err += " dir";
		else
			err += " rc";
		err += " : ";
		err += "  genome_pos: " + to_string(i_genome) + " ";
		err.push_back(genome[i_genome]);
		err += "   seq_pos: " + to_string(i_read) + " ";
		err.push_back(read[i_read]);
		err += "\n";
		for (auto p = md.ext_cigar; *p; ++p)
			err.push_back(*p);
		err += "\n";

		for (auto c : genome)
			err.push_back(c);
		err += "\n";

		for (int i = 0; i < sequence_len; ++i)
			err.push_back(read[i]);
		err += "\n";

		cerr << err;

		return false;
	}

	return true;
}
#endif

#ifdef REFINE_SCORE_SNP
// ************************************************************************************
void CSamGenerator::refine_score_snp(mapping_desc_t& md, uchar_t* id, uint32_t id_len, uchar_t* sequence, uint32_t sequence_len, genome_t direction)
{
	if (!params->enable_var_indel_long && !params->enable_var_indel_short && !params->enable_var_snp)
		return;

	auto ref = reference->GetData();

	vector<uchar_t> genome;
	vector<uchar_t> read;
	vector<array<float, 4>> ref_snp;

	uint32_t begin;
	uint32_t size = sequence_len + 2 * params->max_approx_indel_len;
	genome.resize((size_t)sequence_len + 2 * params->max_approx_indel_len + 2ll);
	ref_snp.resize(genome.size());

	ref_seq_desc->RevTranslate(begin, md.ref_seq_pos, md.ref_seq_id);

	// Genome prefetch
	uint32_t i_genome = 0;
	uchar_t* genome_ptr = ref + (begin >> 1);
	uint32_t text_size_div2 = size / 2;

	if (begin & 1)
		genome[i_genome++] = *genome_ptr++ & 0x0f;
	for (uint32_t i = 0; i < text_size_div2; ++i)
	{
		genome[i_genome++] = *genome_ptr >> 4;
		genome[i_genome++] = *genome_ptr++ & 0x0f;
	}

	// Load SNP data and prepare reference sequence with applied known SNP frequencies
	for (int i = 0; i < ref_snp.size(); ++i)
	{
		for (int j = 0; j < 4; ++j)
			ref_snp[i][j] = 0.0;
		
		if(genome[i] < 4)
			ref_snp[i][genome[i]] = 1.0;
	}

	// Appply SNP data
	auto snp_range = variant_db->FindRangeSNP(begin, begin + ref_snp.size());
	for (auto p = snp_range.first; p != snp_range.second; ++p)
		if(p->symbol < 4 && p->pos - begin < ref_snp.size())
			ref_snp[p->pos - begin][p->symbol] = ((float)p->freq) / 65535.0;

	// Load read sequence
	read.resize(sequence_len);
	if (direction == genome_t::direct)
		for (int i = 0; i < sequence_len; ++i)
		{
			switch (sequence[i])
			{
			case 'A':
				read[sequence_len - 1 - i] = 0; break;
			case 'C':
				read[sequence_len - 1 - i] = 1; break;
			case 'G':
				read[sequence_len - 1 - i] = 2; break;
			case 'T':
				read[sequence_len - 1 - i] = 3; break;
			default:
				read[sequence_len - 1 - i] = 4;
			}
		}
	else
		for (int i = 0; i < sequence_len; ++i)
		{
			switch (sequence[i])
			{
			case 'A':
				read[sequence_len - 1 - i] = 3; break;
			case 'C':
				read[sequence_len - 1 - i] = 2; break;
			case 'G':
				read[sequence_len - 1 - i] = 1; break;
			case 'T':
				read[sequence_len - 1 - i] = 0; break;
			default:
				read[sequence_len - 1 - i] = 4;
			}
		}

	i_genome = 0;
	uint32_t i_read = 0;
	uchar_t* ec = md.ext_cigar;

	double cost_update = 0.0;
	double match_score = params->match_score;
	double mismatch_score = params->mismatch_score;

	while (*ec)
	{
		if (*ec == '.')
		{
			cost_update -= match_score;
			if (read[i_read] < 4)
			{
				cost_update += match_score * ref_snp[i_genome][read[i_read]];
				cost_update += mismatch_score * (1.0 - ref_snp[i_genome][read[i_read]]);
			}
			++i_genome;
			++i_read;
			++ec;
		}
		else if (*ec == 'A' || *ec == 'C' || *ec == 'G' || *ec == 'T' || *ec == 'N')
		{
			cost_update -= mismatch_score;
			if (read[i_read] < 4)
			{
				cost_update += match_score * ref_snp[i_genome][read[i_read]];
				cost_update += mismatch_score * (1.0 - ref_snp[i_genome][read[i_read]]);
			}
			++i_genome;
			++i_read;
			++ec;
		}
		else if (*ec == '$')
		{
			++i_read;
			++ec;
		}
		else if (*ec == '#')
		{
			++ec;
			++i_read;
			++ec;
		}
		else if (*ec == '^')
		{
			++ec;
			++i_genome;
			++ec;
		}
	}

	md.score += cost_update;
}
#endif
// EOF
