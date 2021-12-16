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

#include "mapper.h"
#include <sstream>
#include "../libs/asmlib.h"
#include "../common/utils.h"


// ************************************************************************************
CMapper::CMapper()
{
	// Internal parameters

	// Parameters of mapping
	params.min_read_len = 101;			// minimal length of mapped reads
	params.max_read_len = 101;			// maximal length of mapped reads
	params.max_fastq_rec_length = max_fastq_rec_length;	// 10240
	params.max_cigar_len = 4096;

	objects.running_stats = nullptr;

	objects.serial_processing = nullptr;

	objects.progress = new CProgress();

	total_input_file_size = 0;
}

// ************************************************************************************
CMapper::~CMapper()
{
	delete objects.progress;
}

// ************************************************************************************
bool CMapper::SetParams(const CCmdParams &cmd_params)
{
	params = cmd_params;

	return true;
}

// ************************************************************************************
bool CMapper::SetRunningStats(CRunningStats *_running_stats)
{
	objects.running_stats = _running_stats;

	return true;
}

// ************************************************************************************
bool CMapper::StartMapping()
{
	CStopWatch watch;
	watch.StartTimer();

	if (params.input_file_names.empty() || params.index_name == "" || params.project_name == "")
		return false;

	if (!check_input_files())
		return false;

	adjust_threads(true);
	if (!adjust_bins())
		return false;

	bool res;

	prepare_running_stats(true);

	// Split the reads before mapping
	adjust_memory_splitting();

	reads_splitting();

	prepare_running_stats(false);

	CStopWatch watch_main;
	watch_main.StartTimer();

	// Set the places of segments for substages
	adjust_stage_segments();

	prepare_serial_processing(params.io_calls_in_serial_mode);

	prepare_reference_sa();

	adjust_memory_mapping();

	if (params.mapping_mode == mapping_mode_t::first)
		res = reads_mapping_first_stratum();
	else if (params.mapping_mode == mapping_mode_t::second)
		res = reads_mapping_second_stratum();
	else
		res = reads_mapping_all_strata();

	release_reference_sa();
	release_serial_processing();

	watch_main.StopTimer();
#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_MAIN, watch_main.GetElapsedTime());
#endif
	if (params.verbosity_level > 0)
		cerr << "Main processing time: " << watch_main.GetElapsedTime() << "s\n";

	if (!params.developer_mode && res)
	{
		CStopWatch watch_post;
		watch_post.StartTimer();

		prepare_serial_processing(params.io_calls_in_serial_mode);

		if (params.verbosity_level > 0)
		{
			cerr << "***** Postprocessing *****\n";
			fflush(stderr);
		}

		prepare_reference();
		adjust_memory_postprocessing();
		res = reads_postprocessing(params.mapping_mode);
		release_reference();

		release_serial_processing();

		watch_post.StopTimer();

#ifdef COLLECT_STATS
		objects.running_stats->AddTotals(STAT_TIME_POST, watch_post.GetElapsedTime());
#endif

		if (params.verbosity_level > 0)
			cerr << "Postprocessing time: " << watch_post.GetElapsedTime() << "s\n";
	}

	watch.StopTimer();

#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_TOTAL, watch.GetElapsedTime());
#endif

	return res;
}

// ************************************************************************************
// Perform only a selected stage of mapping - for development purposes only!
bool CMapper::StartMapping(uint32_t stage_major, uint32_t stage_minor, bool sensitive_mode)
{
	if (params.input_file_names.empty() || params.index_name == "" || params.project_name == "")
		return false;

	if (!check_input_files())
		return false;

	adjust_threads(true);
	if (!adjust_bins())
		return false;

	bool res;
	CStopWatch watch;
	watch.StartTimer();
	prepare_running_stats(false);
	prepare_running_stats(true);

	prepare_serial_processing(params.io_calls_in_serial_mode);

	// Set the places of segments for substages
	adjust_stage_segments();

	// Split the reads before mapping
	if (stage_major == 0)
	{
		adjust_memory_splitting();
		res = reads_splitting();
	}
	else
	{
		prepare_reference_sa();
		adjust_memory_mapping();
		res = reads_mapping_single_stage(stage_major, stage_minor, params.mapping_mode, sensitive_mode);
		release_reference_sa();
	}

	watch.StopTimer();
#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_TOTAL, watch.GetElapsedTime());
#endif
	release_serial_processing();

	return res;
}

// ************************************************************************************
bool CMapper::StartMapping(uint32_t stage_major_from, uint32_t stage_minor_from, uint32_t stage_major_to, uint32_t stage_minor_to, bool sensitive_mode)
{
	if (params.input_file_names.empty() || params.index_name == "" || params.project_name == "")
		return false;

	if (!check_input_files())
		return false;

	adjust_threads(true);
	if (!adjust_bins())
		return false;

	bool res = true;
	CStopWatch watch;
	watch.StartTimer();
	prepare_running_stats(false);
	prepare_running_stats(true);

	prepare_serial_processing(params.io_calls_in_serial_mode);

	// Set the places of segments for substages
	adjust_stage_segments();

	// Split the reads before mapping
	if (stage_major_from == 0)
	{
		adjust_memory_splitting();
		res = reads_splitting();
		stage_major_from = 1;
		stage_minor_from = 0;
	}

	if (stage_major_from < stage_major_to || (stage_major_from == stage_major_to && stage_minor_from <= stage_minor_to))
	{
		prepare_reference_sa();
		adjust_memory_mapping();
		res = reads_mapping_stage_range(stage_major_from, stage_minor_from, stage_major_to, stage_minor_to, params.mapping_mode, sensitive_mode);
		release_reference_sa();
	}

	watch.StopTimer();
#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_TOTAL, watch.GetElapsedTime());
#endif
	release_serial_processing();

	return res;
}

// ************************************************************************************
bool CMapper::StartPostProcessing()
{
	if (params.input_file_names.empty() || params.index_name == "" || params.project_name == "")
		return false;

	if (!check_input_files())
		return false;

	adjust_threads(false);
	if (!adjust_bins())
		return false;

	bool res;
	CStopWatch watch;
	watch.StartTimer();
	prepare_running_stats(true);
	prepare_running_stats(false);

	prepare_serial_processing(params.io_calls_in_serial_mode);

	if (params.verbosity_level > 0)
	{
		cerr << "***** Postprocessing *****\n";
		fflush(stderr);
	}

	prepare_reference();
	adjust_memory_postprocessing();
	
	res = reads_postprocessing(params.mapping_mode);
	release_reference();

	watch.StopTimer();
#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_TOTAL, watch.GetElapsedTime());
#endif
	release_serial_processing();

	return res;
}

// ************************************************************************************
// Adjust the number of threads of algorithm parts
bool CMapper::adjust_threads(bool single_thr_readers)
{
	if (params.no_threads == 0)
		params.no_threads = thread::hardware_concurrency();
	if (params.no_threads == 0)		// due to failure of thread::hardware_concurrency()
		params.no_threads = 1;

	// Splitting stage
	params.no_thr_fastq_reader = MAX(1, params.no_threads / 2);
	uint32_t no_input_files = 0;

	if (single_thr_readers)
	{
		for (auto &x : params.input_file_names)
			if (x.second.empty())
				no_input_files += 1;
			else
				no_input_files += 2;
	}
	else
		no_input_files = (uint32_t) (params.input_file_names.size());

	if (!params.input_file_names.empty() && no_input_files < params.no_thr_fastq_reader)
		params.no_thr_fastq_reader = no_input_files;
	params.no_thr_splitter = MAX(params.no_threads - params.no_thr_fastq_reader + 1, 1);

	// Mapping stage
	params.no_thr_mapping_cores = params.no_threads;

	// SAM production stage
	params.no_thr_sam_generators = MAX(params.no_threads - params.no_threads / 24, 1);

	return true;
}

// ************************************************************************************
// Adjust the size of memory for queues in splitting stage
bool CMapper::adjust_memory_splitting()
{
	int64_t mem_available = params.max_total_memory;

	// Memory for FASTQ blocks
	params.block_overhead_size = 1 << 16;
	if (params.paired_reads)
		params.block_size = 1 << 25;
	else
		params.block_size = 1 << 25;

	for (params.no_fastq_blocks_per_thread = 12; params.no_fastq_blocks_per_thread > 3; --params.no_fastq_blocks_per_thread)
		if ((params.block_size + params.block_overhead_size) * params.no_thr_fastq_reader * params.no_fastq_blocks_per_thread < mem_available * 0.25)
			break;

	mem_available -= (params.block_size + params.block_overhead_size) * params.no_thr_fastq_reader * params.no_fastq_blocks_per_thread;
	if (mem_available <= 0)
	{
		cerr << "Not enough memory\n";
		exit(1);
	}

	// Memory for reads
	mem_available -= params.read_len * params.no_thr_splitter * 4;

	// Memory for bin writer
	params.max_bin_part_size = 1 << 27;
	if (params.max_bin_part_size > 0.05 * mem_available)
		params.max_bin_part_size = (uint64_t)(0.05 * mem_available);
	mem_available -= params.max_bin_part_size;

	// Memory for bins in splitting stage
	if (params.ssd_mode)
		params.bin_size = 1 << 15;
	else
		params.bin_size = 1 << 16;

	for (; params.bin_size > 128; params.bin_size = (uint64_t)(params.bin_size * 0.75))
	{
		params.no_parts = (mem_available / 2) / params.bin_size;
		if (params.no_parts > 2ull * params.no_bins * params.no_thr_splitter + params.no_parts / 16)
			break;
	}

	params.max_bin_part_size = 1 << 27;
	params.no_parts = (mem_available / params.bin_size);
	params.min_no_free_parts = params.no_parts / 16;

	if (params.verbosity_level > 1)
	{
		cerr << "Max total memory             : " << FormatInt(params.max_total_memory) << "\n";
		cerr << "Memory settings fo splitting stage:\n";
		cerr << "  no. FASTQ blocks per thread: " << FormatInt(params.no_fastq_blocks_per_thread) << "\n";
		cerr << "  mem. per FASTQ blocks      : " << FormatInt((params.block_size + params.block_overhead_size) * params.no_thr_fastq_reader * params.no_fastq_blocks_per_thread) << "\n";
		cerr << "  max bin part size          : " << FormatInt(params.max_bin_part_size) << "\n";
		cerr << "  bin size                   : " << FormatInt(params.bin_size) << "\n";
		cerr << "  no bin parts               : " << FormatInt(params.no_parts) << "\n";
		cerr << "  min. no. free bin parts    : " << FormatInt(params.min_no_free_parts) << "\n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
// Adjust the size of memory for queues
bool CMapper::adjust_memory_mapping()
{
	int64_t mem_available = params.max_total_memory;

	// Subtract the amount of memory for reference sequence
	mem_available -= objects.reference->GetSize();

	// Subtract estimated amount of memory for SA parts
	uint64_t no_parts = params.no_bins * (1ull << (params.sa_prefix_overhead * 2));
	uint64_t mem_for_SA_parts = objects.reference->GetSize() / no_parts * params.no_thr_mapping_cores; // no. of contemporarily used parts 
	mem_for_SA_parts *= 2;				// count both SA dir and SA r.c.
	mem_for_SA_parts *= 4;				// 4B for each element
	mem_for_SA_parts *= 2;				// double the space for safety
	mem_available -= mem_for_SA_parts;

	if (mem_available < (1ll << 28))
	{
		cerr << "Not enough memory\n";
		exit(1);
	}

	// Memory for input bins
	params.max_mapping_memory = mem_available / 2;
	mem_available -= params.max_mapping_memory;
	mem_available -= params.read_len * params.no_thr_mapping_cores;

	// Memory for bin writer
	params.max_bin_part_size = 1 << 27;
	if (params.max_bin_part_size > 0.05 * mem_available)
		params.max_bin_part_size = (uint64_t)(0.05 * mem_available);
	mem_available -= params.max_bin_part_size;

	// Memory for results
	for (params.res_group_size = 1 << 15; params.res_group_size > 128; params.res_group_size = (uint64_t)(params.res_group_size * 0.75))
	{
		params.no_res_parts = (uint32_t)((mem_available / 2) / params.res_group_size);
		if (params.no_res_parts > 2 * params.no_res_groups * params.no_thr_mapping_cores + params.no_res_parts / 16)
			break;
	}

	params.max_group_part_size = 1 << 27;
	params.no_res_parts = (uint32_t)((mem_available / 2) / params.res_group_size);
	params.min_no_free_group_parts = params.no_res_parts / 16;

	// Memory for bins in splitting stage
	if (params.ssd_mode)
		params.bin_size = 1 << 15;
	else
		params.bin_size = 1 << 16;

	// Memory for output bins
	for (; params.bin_size > 128; params.bin_size = (uint64_t)(params.bin_size * 0.75))
	{
		params.no_parts = (mem_available / 2) / params.bin_size;
		if (params.no_parts > 2ull * params.no_bins * params.no_thr_mapping_cores + params.no_parts / 16)
			break;
	}

	params.max_bin_part_size = 1 << 27;
	params.no_parts = (mem_available / 2) / params.bin_size;
	params.min_no_free_parts = params.no_parts / 16;

	if (params.verbosity_level > 1)
	{
		cerr << "Memory settings fo mapping stage:\n";
		cerr << "  max mapping memory       : " << FormatInt(params.max_mapping_memory) << "\n";
		cerr << "  no. result parts         : " << FormatInt(params.no_res_parts) << "\n";
		cerr << "  result group size        : " << FormatInt(params.res_group_size) << "\n";
		cerr << "  min. no. free group parts: " << FormatInt(params.min_no_free_group_parts) << "\n";
		cerr << "  max bin part size        : " << FormatInt(params.max_bin_part_size) << "\n";
		cerr << "  bin part size            : " << FormatInt(params.bin_size) << "\n";
		cerr << "  no. bin parts            : " << FormatInt(params.no_parts) << "\n";
		cerr << "  min. no. free bin parts  : " << FormatInt(params.min_no_free_parts) << "\n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
// Adjust memory for post processing stage
bool CMapper::adjust_memory_postprocessing()
{
	int64_t mem_available = params.max_total_memory;

	// Memory for reference
	mem_available -= ref_seq_desc.GetSize(params.index_name);

	// Memory for reads
	mem_available -= params.read_len * params.no_thr_splitter * 4;

	mem_available -= (int64_t) params.max_fastq_rec_length * params.no_thr_sam_generators * 3 * 4;

	// Memory for CIGAR and MD
	mem_available -= 6ll * params.max_cigar_len * params.no_thr_sam_generators * 4 * params.max_no_mappings;

	// Memory for SAM parts
	params.sam_buffer_memory = 128ull << 20;
	mem_available -= 2ll * params.sam_buffer_memory;

	if (params.store_BAM)
		params.sam_part_size = 64ull << 10;
	else
		params.sam_part_size = 64ull << 10;
//		params.sam_part_size = 16 << 20;

	int64_t mem_sam_parts;

	if (params.gzipped_SAM_level == 0)
		mem_sam_parts = params.no_thr_sam_generators * params.sam_part_size * 2;
	else
		mem_sam_parts = params.no_thr_sam_generators * params.sam_part_size * 3;

	mem_available -= mem_sam_parts;

	if (mem_available <= 0)
	{
		cerr << "Not enough memory\n";
		exit(1);
	}

	// Memory for FASTQ blocks - must be of the same size as in splitting stage!!!
	params.block_overhead_size = 1 << 16;
	if (params.paired_reads)
		params.block_size = 1 << 25;
	else
		params.block_size = 1 << 25;

	params.no_fastq_blocks_per_thread = 1ull << params.id_bits_subgroup;

/*	for (params.no_fastq_blocks_per_thread = (1 << params.id_bits_subgroup) * 4;
		params.no_fastq_blocks_per_thread > (1 << params.id_bits_subgroup) * 2.3;
		--params.no_fastq_blocks_per_thread)
	{
		if (params.verbosity_level > 1)
		{
			cerr << "No. fastq blocks " << params.no_fastq_blocks_per_thread << "      " << params.id_bits_subgroup << "\n";
			cerr << (params.block_size + params.block_overhead_size) * params.no_fastq_blocks_per_thread * (params.no_thr_fastq_reader + params.no_thr_sam_generators) << "\n";
			cerr << params.block_size << "  " << params.block_overhead_size << "  " << params.no_fastq_blocks_per_thread << "  " << params.no_thr_fastq_reader << "  " << params.no_thr_sam_generators << "\n";
			cerr << mem_available * 0.15 << "\n";

			fflush(stderr);
		}
//		if ((params.block_size + params.block_overhead_size) * params.no_fastq_blocks_per_thread * (params.no_thr_fastq_reader + params.no_thr_sam_generators) < (1ull << 30))
		if ((params.block_size + params.block_overhead_size) * params.no_fastq_blocks_per_thread * (params.no_thr_fastq_reader + params.no_thr_sam_generators) < (1ull << 30))
			break;
	}

	if (params.verbosity_level > 1)
	{
		cerr << "No. fastq blocks " << params.no_fastq_blocks_per_thread << "   " << (1 << params.id_bits_subgroup) << "\n";
		fflush(stderr);
	}
	*/
/*	mem_available -= (params.block_size + params.block_overhead_size) * (params.no_thr_fastq_reader + params.no_thr_sam_generators) * params.no_fastq_blocks_per_thread;
	if (mem_available <= 0)
	{
		cerr << "Not enough memory\n";
		exit(1);
	}*/

	// Memory for result groups
//	if (mem_available < ((int64_t) params.no_threads) * (1024ll << 20))
	mem_available = ((int64_t) params.no_threads) * (1024ll << 20);
	params.max_res_dev_memory = mem_available;

	if (params.verbosity_level > 1)
	{
		cerr << "Max total memory             : " << FormatInt(params.max_total_memory) << "\n";
		cerr << "Memory settings for postprocessing stage:\n";
		cerr << "  no. FASTQ blocks per thread: " << FormatInt(params.no_fastq_blocks_per_thread) << "\n";
		cerr << "  mem. per FASTQ blocks      : " << FormatInt(((uint64_t) params.block_size + params.block_overhead_size) * 
			((uint64_t) params.no_thr_fastq_reader + params.no_thr_sam_generators) * params.no_fastq_blocks_per_thread) << "\n";
		cerr << "  mem for groups delivery    : " << FormatInt(params.max_res_dev_memory) << "\n";
		cerr << "  mem for SAM parts          : " << FormatInt(mem_sam_parts) << "\n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
// Read the suffix array of the reference sequence and adjust the bins to make the distribution
// of sizes of bins as flat as possible
bool CMapper::adjust_bins()
{
#ifdef _DEBUG 
	uint32_t max_prefix_len = 5;		// maximal length of prefix
#else
	uint32_t max_prefix_len = 10;		// maximal length of prefix
	//uint32_t max_prefix_len = 5;		// maximal length of prefix
#endif

	shared_ptr<CSuffixArray> sa(new CSuffixArray());
	if (!sa->SetIndexName(params.index_name, genome_t::direct))
	{
		cerr << "Error: Cannot open index files\n";
		return false;
	}

	uint32_t *lut_data;
	uint32_t lut_prefix_len;

	sa->GetShortLUT(lut_data, lut_prefix_len);

	CBinPrefixes bpx(max_prefix_len, params.verbosity_level);
	bpx.SetLUT(lut_data, lut_prefix_len);
	bpx.ReduceMap(params.no_bins);

	bpx.GetMaps(params.prefix_map, params.prefix_file_names);

	return true;
}

// ************************************************************************************
// Adjust splitting of reads for the stages
bool CMapper::adjust_stage_segments()
{
	uint32_t max_errors = params.max_no_errors;

	params.stage_segments.clear();
	params.stage_segments.resize(max_errors + 2ull);

	for (uint32_t i = 0; i <= max_errors; ++i)
		for (uint32_t j = 0; j <= i; ++j)
			params.stage_segments[i].push_back(make_pair((params.min_read_len*j + i) / (i + 1), (params.min_read_len*(j + 1) + i) / (i + 1)));

	// ``Sentinel''
	params.stage_segments[max_errors + 1ull].push_back(make_pair(0, params.min_read_len));

	return true;
}

// ************************************************************************************
uint32_t CMapper::get_next_stage(uint32_t c_stage, uint32_t m_stage)
{
	uint32_t r_stage;

	if (c_stage == m_stage)
		return m_stage + 1;

	if (m_stage <= 3)
		r_stage = c_stage + 1;
	else if (m_stage <= 5)
	{
		if (c_stage < 2)
			r_stage = c_stage + 1;
		else
			r_stage = m_stage;
	}
	else if (m_stage == 6)
	{
		if (c_stage < 3)
			r_stage = c_stage + 1;
		else
			r_stage = m_stage;
	}
	else if (m_stage == 7)
	{
		if (c_stage < 4)
			r_stage = c_stage + 1;
		else
			r_stage = m_stage;
	}
	else
	{
		if (c_stage < m_stage - 3)
			r_stage = c_stage + 1;
		else
			r_stage = m_stage;
	}

	return r_stage;
}

// ************************************************************************************
// Adjust maximal and minimal read length (if not given explicitely in mapper execution)
// Also adjust maximal no. of errors as a fraction of min. read length
bool CMapper::adjust_max_and_min_read_len(vector<uint64_t> &hist_read_len)
{
	auto p = find_if(hist_read_len.begin(), hist_read_len.end(), [](uint64_t x) {return x != 0; });
	auto q = find_if(hist_read_len.rbegin(), hist_read_len.rend(), [](uint64_t x) {return x != 0; });

	if (p == hist_read_len.end())		// No reads
		return false;

	int min_len = (int) (p - hist_read_len.begin());
	int max_len = (int) (hist_read_len.rend() - q - 1);

	if (min_len < (int) shortest_min_read_len)
		min_len = (int) shortest_min_read_len;
	if (min_len > max_len)
		max_len = min_len;

	params.max_read_len = max_len;
	params.min_read_len = max_len;

	params.max_cigar_len = params.max_read_len * 5;

	if (min_len < max_len)
	{
		// Allow min read len to be 10% smaller than max read len
		params.min_read_len = MAX(min_len, (int)(max_len * 0.9));
	}

	if (params.max_no_errors == 0)
	{
		// Adjust max. no of errors according to read lengths
		if (params.max_read_len < 100)
			params.max_no_errors = 3;
		else if (params.max_read_len < 127)
			params.max_no_errors = 4;
		else if (params.max_read_len < 153)
			params.max_no_errors = 5;
		else if (params.max_read_len < 203)
			params.max_no_errors = 6;
		else if (params.max_read_len < 253)
			params.max_no_errors = 7;
		else
			params.max_no_errors = 8;

		params.max_approx_indel_mismatches = (uint32_t)(params.max_no_errors * params.sensitivity_factor);
	}

	double max_frac_errors = (double)params.max_no_errors / params.min_read_len;

	if(max_frac_errors > largest_frac_errors)
		params.max_no_errors = (uint32_t)(params.min_read_len * largest_frac_errors);

	return true;
}

// ************************************************************************************
// Prepare reference and SA
bool CMapper::prepare_reference_sa()
{
	objects.reference = new CReference();
	objects.sa_dir = new CSuffixArray(params.sa_in_ram);
	objects.sa_rc = new CSuffixArray(params.sa_in_ram);

#ifdef ENABLE_VCF_VARIANTS
	objects.variant_db = new CVariantDB();
#endif

	bool res = true;

	// Load reference
	if (params.verbosity_level > 0)
	{
		cerr << "** Loading reference and index **\n";
		fflush(stderr);
	}
	res &= objects.reference->SetIndexName(params.index_name);

	// Prepare for SA load
	if (params.verbosity_level > 1)
	{
		cerr << "Preparing suffix arrays\n";
		fflush(stderr);
	}
	res &= objects.sa_dir->SetIndexName(params.index_name, genome_t::direct);
	res &= objects.sa_rc->SetIndexName(params.index_name, genome_t::rev_comp);

#ifdef ENABLE_VCF_VARIANTS
	// Load variants (if required)
	if (params.enable_var_indel_long || params.enable_var_indel_short || params.enable_var_snp)
	{
		if (params.verbosity_level > 0)
		{
			cerr << "** Loading variants database **\n";
			fflush(stderr);
		}

		shared_ptr<CMapperFile> vf_vcf(new CMapperFile(EXT_VCF, MARKER_VCF));
		shared_ptr<CMapperFile> vf_ref_snp(new CMapperFile(EXT_REF_SNP, MARKER_REF_SNP));

		if(vf_vcf->OpenRead(params.index_name) && vf_ref_snp->OpenRead(params.index_name))
			objects.variant_db->Deserialize(vf_vcf, vf_ref_snp);
		else
			cerr << "Variant file index does not exist\n";
	}
#endif

	// Prepare ID store
	CIDStore *id_store = new CIDStore(params.id_bits_total, params.id_bits_subgroup, params.id_bits_local, params.verbosity_level);
	shared_ptr<CMapperFile> mf(new CMapperFile(EXT_ID_STORE, MARKER_ID_STORE));
	if (!mf->OpenRead(params.temp_prefix + PRJ_CFG_FILE))
	{
		cerr << "Cannot open project configuration file: " << params.temp_prefix + PRJ_CFG_FILE << "\n";
		exit(1);
	}

	id_store->Load(mf);
	params.no_res_groups = id_store->GetNumGroups();
	mf->Close();

	delete id_store;

	return res;
}

// ************************************************************************************
// Prepare reference only
bool CMapper::prepare_reference()
{
	objects.reference = new CReference();
#ifdef ENABLE_VCF_VARIANTS
	objects.variant_db = new CVariantDB();
#endif

	bool res = true;

	if (params.verbosity_level > 0)
	{
		cerr << "** Loading reference **\n";
		fflush(stderr);
	}

	res = objects.reference->SetIndexName(params.index_name);

#ifdef ENABLE_VCF_VARIANTS
	// Load variants (if required)
	if (params.enable_var_indel_long || params.enable_var_indel_short || params.enable_var_snp)
	{
		if (params.verbosity_level > 0)
		{
			cerr << "** Loading variants database **\n";
			fflush(stderr);
		}

		shared_ptr<CMapperFile> vf_vcf(new CMapperFile(EXT_VCF, MARKER_VCF));
		shared_ptr<CMapperFile> vf_ref_snp(new CMapperFile(EXT_REF_SNP, MARKER_REF_SNP));

		if (vf_vcf->OpenRead(params.index_name) && vf_ref_snp->OpenRead(params.index_name))
			objects.variant_db->Deserialize(vf_vcf, vf_ref_snp);
		else
			cerr << "Variant file index does not exist\n";
	}
#endif

	return res;
}

// ************************************************************************************
// Prepare reference and SA
bool CMapper::release_reference_sa()
{
	delete objects.reference;
	delete objects.sa_dir;
	delete objects.sa_rc;

	objects.reference = nullptr;
	objects.sa_dir = nullptr;
	objects.sa_rc = nullptr;

#ifdef ENABLE_VCF_VARIANTS
	delete objects.variant_db;
	objects.variant_db = nullptr;
#endif

	return true;
}

// ************************************************************************************
// Prepare reference only
bool CMapper::release_reference()
{
	delete objects.reference;
	objects.reference = nullptr;

#ifdef ENABLE_VCF_VARIANTS
	delete objects.variant_db;
	objects.variant_db = nullptr;
#endif

	return true;
}

// ************************************************************************************
// Prepare object for serial processing
bool CMapper::prepare_serial_processing(bool call_in_serial_mode)
{
	objects.serial_processing = new CSerialProcessing(call_in_serial_mode);

	return true;
}

// ************************************************************************************
// Release object for serial processing
bool CMapper::release_serial_processing()
{
	if (objects.serial_processing)
		delete objects.serial_processing;

	return true;
}

// ************************************************************************************
// Prepare running stats
bool CMapper::prepare_running_stats(bool before_preprocessing)
{
	if (!objects.running_stats)
		return false;

#ifdef COLLECT_STATS
	if (before_preprocessing)
	{
		objects.running_stats->Register(1000001, "Fastq Readers, send bytes: ", running_stats_t::lists);

		objects.running_stats->Register(STAT_SORTING_TOTAL, "Reads sorting (total)", running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_TOTAL, "Total computation time: total", running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_SPLIT, "Total computation time: reads split", running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_MAIN, "Total computation time: main processing", running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_POST, "Total computation time: post processing", running_stats_t::totals);

		objects.running_stats->Register(STAT_DIR_EXACT_MATCH, "No. of reads with DIR exact match: ", running_stats_t::totals);
		objects.running_stats->Register(STAT_RC_EXACT_MATCH, "No. of reads with RC exact match: ", running_stats_t::totals);
		objects.running_stats->Register(STAT_DIR_AND_RC_EXACT_MATCH, "No. of reads with DIR_AND_RC exact match: ", running_stats_t::totals);
		objects.running_stats->Register(STAT_DIR_OR_RC_EXACT_MATCH, "No. of reads with DIR_OR_RC exact match: ", running_stats_t::totals);

		objects.running_stats->Register(STAT_PUSH_RESULTS, "Push results ", running_stats_t::lists);
		objects.running_stats->Register(STAT_PUSH_RESULTS_UNQ, "Push results unq. ", running_stats_t::lists);
		objects.running_stats->Register(STAT_PUSH_RESULTS_SEND_BYTES, "Push results - send bytes ", running_stats_t::averages);
		objects.running_stats->Register(STAT_PUSH_RESULTS_ADDED_BYTES, "Push results - added bytes ", running_stats_t::averages);

		objects.running_stats->Register(STAT_MAPPED_SE_SINGLE, "No. of mapped SE reads", running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_SE_NONE, "No. of unmapped SE reads", running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_PE_PAIR_UNQ, "No. of mapped PE reads - pair (1 mapping)", running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_PE_PAIR_MORE, "No. of mapped PE reads - pair (>1 mapping)", running_stats_t::totals);

		objects.running_stats->Register(STAT_MAPPED_PE_ERRORS, "No. of mapped PE reads - paired errors (single reads)", running_stats_t::totals);

		objects.running_stats->Register(STAT_MAPPED_PE_INDEPENDENT_BOTH, "No. of mapped PE reads - unpaired: both mapped", running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_PE_INDEPENDENT_SINGLE, "No. of mapped PE reads - unpaired: one mapped", running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_PE_INDEPENDENT_ERRORS, "No. of mapped PE reads - unpaired errors", running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_PE_INDEPENDENT_NONE, "No. of unmapped PE reads", running_stats_t::totals);

		objects.running_stats->Register(STAT_SHORT_INDEL_REFINEMENTS_LEV, "No. of refined short indels from Lev", running_stats_t::totals);
		objects.running_stats->Register(STAT_SHORT_INDEL_REFINEMENTS_MAPPING, "No. of refined short indels from mapping", running_stats_t::totals);

		for (int i = 0; i < (int) MatchingMethod::size(); ++i) {
			ostringstream oss;
			MatchingMethod mm((MatchingMethod::Enumeration)i);
			oss << "No. of mapped PE reads - " << mm.toString();
			objects.running_stats->Register(STAT_MAPPED_PE_METHOD + i, oss.str(), running_stats_t::totals);
		}

		for (int i = 0; i < 51; ++i) {
			ostringstream oss;
			oss << "PE reads histogram - " << setfill('0') << setw(2) << i << " errors";
			objects.running_stats->Register(STAT_MAPPED_PE_HISTO + i, oss.str(), running_stats_t::totals);
		}

		for (int i = 0; i < 51; ++i) {
			ostringstream oss;
			oss << setfill('0') << setw(3) << "PE reads clipping histogram - [" << i * 5 << "," << (i + 1) * 5 << ")";
			objects.running_stats->Register(STAT_MAPPED_PE_HISTO_CLIPPED + i, oss.str(), running_stats_t::totals);
		}
	}
	else
	{
		objects.running_stats->Register(1000002, "SAM Generators, read bytes: ", running_stats_t::lists);

		uint32_t max_errors = params.max_no_errors;
		for (uint32_t i = 0; i <= max_errors; ++i)
			for (uint32_t j = 0; j <= i; ++j)
			{
				objects.running_stats->Register(STAT_SORTING_BASE + (i << 5) + j, "Reads sorting: stage " + StageDesc(i, j), running_stats_t::averages);
				objects.running_stats->Register(STAT_TIME_BASE + (i << 5) + j, "Computation time: stage " + StageDesc(i, j), running_stats_t::lists);

				objects.running_stats->Register(STAT_MISMATCHES_TO_DIR_BASE + (i << 5) + j, "No. of reads matching DIR gen at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_MISMATCHES_TO_RC_BASE + (i << 5) + j, "No. of reads matching RC gen at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_SELECTED_READS_BASE + (i << 5) + j, "No. of SELECTED reads at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_ACCEPTED + (i << 5) + j, "No. of TESTED reads (cf+) at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_DISCARDED + (i << 5) + j, "No. of DISCARDED reads (CF-) at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_ACCEPTED_IN_MALICIOUS_GROUPS + (i << 5) + j, "No. of TESTED reads (cf+) IN MALICIOUS_GROUPS at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_DISCARDED_IN_MALICIOUS_GROUPS + (i << 5) + j, "No. of DISCARDED reads (CF-)IN MALICIOUS_GROUPS at the " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_LEV_POSITIVE + (i << 5) + j, "No. of positive tests (LEV+) made by Lev_diag2 " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_INDEL_MAPPING_FOUND + (i << 5) + j, "No. of found indels during mapping " + StageDesc(i, j), running_stats_t::lists);
				objects.running_stats->Register(STAT_INDEL_MAPPING_TESTS + (i << 5) + j, "No. of examined indels during mapping " + StageDesc(i, j), running_stats_t::lists);
				//objects.running_stats->Register(STAT_LEV_NEGATIVE + (i << 5) + j, "No. of negative tests (LEV-) made by Lev_diag2 " + StageDesc(i, j), running_stats_t::lists);
			}

		if (params.sensitive_mode)
			for (uint32_t j = 0; j <= max_errors; ++j)
			{
				objects.running_stats->Register(STAT_SORTING_BASE + (max_errors << 5) + j + (1 << 10), "Reads sorting: stage " + StageDesc(max_errors, j, true), running_stats_t::averages);
				objects.running_stats->Register(STAT_TIME_BASE + (max_errors << 5) + j + (1 << 10), "Computation time: stage " + StageDesc(max_errors, j, true), running_stats_t::lists);

				objects.running_stats->Register(STAT_MISMATCHES_TO_DIR_BASE + (max_errors << 5) + j + (1 << 10), "No. of reads matching DIR gen at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_MISMATCHES_TO_RC_BASE + (max_errors << 5) + j + (1 << 10), "No. of reads matching RC gen at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_SELECTED_READS_BASE + (max_errors << 5) + j + (1 << 10), "No. of SELECTED reads at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_ACCEPTED + (max_errors << 5) + j + (1 << 10), "No. of TESTED reads (cf+) at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_DISCARDED + (max_errors << 5) + j + (1 << 10), "No. of DISCARDED reads (CF-) at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_ACCEPTED_IN_MALICIOUS_GROUPS + (max_errors << 5) + j + (1 << 10), "No. of TESTED reads (cf+) IN MALICIOUS_GROUPS at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_CF_DISCARDED_IN_MALICIOUS_GROUPS + (max_errors << 5) + j + (1 << 10), "No. of DISCARDED reads (CF-)IN MALICIOUS_GROUPS at the " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_LEV_POSITIVE + (max_errors << 5) + j + (1 << 10), "No. of positive tests (LEV+) made by Lev_diag2 " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_INDEL_MAPPING_FOUND + (max_errors << 5) + j + (1 << 10), "No. of found indels during mapping " + StageDesc(max_errors, j, true), running_stats_t::lists);
				objects.running_stats->Register(STAT_INDEL_MAPPING_TESTS + (max_errors << 5) + j + (1 << 10), "No. of examined indels during mapping " + StageDesc(max_errors, j, true), running_stats_t::lists);
				//objects.running_stats->Register(STAT_LEV_NEGATIVE + (max_errors << 5) + j + (1 << 10), "No. of negative tests (LEV-) made by Lev_diag2 " + StageDesc(max_errors, j, true), running_stats_t::lists);
			}
	}

	for (uint32_t j = 0; j <= STAT_MAX_LEN_INDELS; ++j)
	{
		objects.running_stats->Register(STAT_MAPPED_PE_HISTO_VAR_INS_LONG + j, "Variant long indel " + to_string(j), running_stats_t::totals);
		objects.running_stats->Register(STAT_MAPPED_PE_HISTO_VAR_DEL_LONG + j, "Variant long indel " + to_string(j), running_stats_t::totals);
	}

	for (uint32_t j = 0; j <= STAT_MAX_READ_LEN; ++j)
	{
		objects.running_stats->Register(STAT_READS_LEN + j, "Read len " + Int2StringFilled(j, 3), running_stats_t::totals);
		objects.running_stats->Register(STAT_READS_LEN_WO_NS + j, "Read len w/o Ns " + Int2StringFilled(j, 3), running_stats_t::totals);
		objects.running_stats->Register(STAT_READS_NS + j, "Read with Ns " + Int2StringFilled(j, 3), running_stats_t::totals);
	}
#endif

	return true;
}

// ************************************************************************************
// Read FASTQ files and split them into bins
bool CMapper::reads_splitting()
{
#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_TIME_THR_SPLITTER, "Thread time: reads splitting     ", running_stats_t::lists);
	objects.running_stats->Register(STAT_TIME_THR_FASTQ_READER, "Thread time: FASTQ reader        ", running_stats_t::lists);
	objects.running_stats->Register(STAT_TIME_THR_BIN_WRITER_BASE, "Thread time: bin writer (split)  ", running_stats_t::totals);
#endif

	objects.progress->Init(total_input_file_size, true);
	objects.progress->SetComment("");

	CStopWatch watch;
	watch.StartTimer();

	if (params.verbosity_level > 0)
	{
		cerr << "***** Preprocessing of reads *****\n";
		fflush(stderr);
	}

	// Prepare pool memory allocators
	objects.mp_fastq_blocks = new CMemoryPool<uchar_t>((params.block_size + params.block_overhead_size) * params.no_thr_fastq_reader * params.no_fastq_blocks_per_thread,
		params.block_size + params.block_overhead_size, "mp_fastq_blocks");
	objects.mp_bins_write = new CMemoryPool<uchar_t>(params.no_parts * params.bin_size, params.bin_size, "mp_bins_write");
	objects.mp_reads = new CMemoryPool<uchar_t>(params.read_len * params.no_thr_splitter * 4, params.read_len, "mp_reads");

	// Prepare queues
	objects.q_file_names = new CRegisteringQueue<file_name_no_t>(1);
	objects.q_blocks = new CRegisteringQueue<fastq_block_t>(params.no_thr_fastq_reader);
	objects.q_bins_write = new CRegisteringQueue<reads_bin_t>(params.no_thr_splitter);

	// Prepare ID store
	CIDStore *id_store = new CIDStore(params.id_bits_total, params.id_bits_subgroup, params.id_bits_local, params.verbosity_level);
	shared_ptr<CMapperFile> mf(new CMapperFile(EXT_ID_STORE, MARKER_ID_STORE));

	// Load input file names to the queue
	uint32_t file_no = 0;
	for (auto &p : params.input_file_names)
	{
		objects.q_file_names->Push(file_name_no_t(file_no, p.first, p.second));
		file_no += 2;
	}
	objects.q_file_names->MarkCompleted();

	// Create bin writer thread
	CBinsWriter *bw = new CBinsWriter(&params, &objects, "stage_0_0_", 0, false);
	thread *thr_bw = new thread(ref(*bw));

	// Create read splitter threads
	uint32_t max_errors = params.max_no_errors;
	vector<CWReadsSplitter*> wrs(params.no_thr_splitter);
	vector<thread*> thr_rs(params.no_thr_splitter);
	for (uint32_t i = 0; i < params.no_thr_splitter; ++i)
	{
		wrs[i] = new CWReadsSplitter(&params, &objects, false, params.mapping_mode != mapping_mode_t::all ? 1 : max_errors);
		thr_rs[i] = new thread(ref(*wrs[i]));
	}

	// Create FASTQ reader threads
	vector<CPreFastqReader*> fqr(params.no_thr_fastq_reader);
	vector<thread*> thr_fqr(params.no_thr_fastq_reader);
	for (uint32_t i = 0; i < params.no_thr_fastq_reader; ++i)
	{
		fqr[i] = new CPreFastqReader(false, &params, &objects, id_store);
		thr_fqr[i] = new thread(ref(*fqr[i]));
	}

	// Wait for completing the stage
	for (auto &p : thr_fqr)
		p->join();

	if (params.verbosity_level > 0)
	{
		cerr << "\nCompleting the preprocessing (could take a minute or so)\n";
		fflush(stderr);
	}

	for (auto &p : thr_rs)
		p->join();
	thr_bw->join();

	mf->OpenWrite(params.temp_prefix + PRJ_CFG_FILE);
	id_store->Save(mf);

	// Destruct objects
	delete thr_bw;
	delete bw;

	vector<uint64_t> hist_read_len, thr_hist_read_len;

	hist_read_len.resize(params.max_fastq_rec_length, 0);

	for (uint32_t i = 0; i < params.no_thr_splitter; ++i)
	{
		delete thr_rs[i];

		wrs[i]->GetHistReadLen(thr_hist_read_len);

		for (uint32_t j = 0; j < params.max_fastq_rec_length; ++j)
			hist_read_len[j] += thr_hist_read_len[j];

		delete wrs[i];
	}
	thr_rs.clear();	thr_rs.shrink_to_fit();
	wrs.clear();	wrs.shrink_to_fit();

	adjust_max_and_min_read_len(hist_read_len);

	for (uint32_t i = 0; i < params.no_thr_fastq_reader; ++i)
	{
		delete thr_fqr[i];
		delete fqr[i];
	}
	thr_fqr.clear();	thr_fqr.shrink_to_fit();
	fqr.clear();		fqr.shrink_to_fit();

	delete objects.mp_fastq_blocks;
	delete objects.mp_bins_write;
	delete objects.mp_reads;
	delete objects.q_file_names;
	delete objects.q_blocks;
	delete objects.q_bins_write;

	params.no_res_groups = id_store->GetNumGroups();
	delete id_store;

	watch.StopTimer();
#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_SPLIT, watch.GetElapsedTime());
#endif

	if (params.verbosity_level > 0)
		cerr << "Preprocessing time: " << watch.GetElapsedTime() << "s\n";

	return true;
}

// ************************************************************************************
// Read bins and perform mapping (first stratum)
bool CMapper::reads_mapping_first_stratum()
{
	CStopWatch watch;

	string prefix_name = "main";

	if (params.verbosity_level > 0)
	{
		cerr << "***** Reads mapping *****\n";
		fflush(stderr);
	}

	uint32_t max_stage = params.max_no_errors;

#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_SORTING_TOTAL, "Reads sorting (total)", running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_RES_WRITER, "Thread time: results writer      ", running_stats_t::totals);
#endif
	// Prepare stage-substage descriptions
	vector<tuple<uint32_t, uint32_t, bool>> stages_to_do;
	for (uint32_t stage_major = 0; stage_major <= max_stage;)
	{
		for (uint32_t stage_minor = 0; stage_minor <= stage_major; ++stage_minor)
			stages_to_do.push_back(make_tuple(stage_major, stage_minor, false));

		stage_major = get_next_stage(stage_major, max_stage);
	}

	if (params.sensitive_mode)
		for (uint32_t stage_minor = 0; stage_minor <= max_stage; ++stage_minor)
			stages_to_do.push_back(make_tuple(max_stage, stage_minor, true));

	// Mapping results queue, thread, and memory pool
	uint32_t tot_substages = (uint32_t)stages_to_do.size() - 1;

	objects.progress->Init((int64_t) tot_substages * params.no_bins, true);

	objects.mp_reads = new CMemoryPool<uchar_t>(params.read_len * params.no_thr_mapping_cores, params.read_len, "mp_reads");
	objects.q_map_res = new CRegisteringQueue<res_group_t>(params.no_thr_mapping_cores * tot_substages);
	objects.mp_map_res = new CMemoryPool<uchar_t>(params.no_res_parts * params.res_group_size,
		params.res_group_size, "mp_map_res");
	CResultGroupsWriter *rgw = new CResultGroupsWriter(&params, &objects, prefix_name);
	thread *thr_rgw = new thread(ref(*rgw));

	for (auto p = stages_to_do.begin() + 1; p != stages_to_do.end(); ++p)
	{
		uint32_t stage_major = get<0>(*p);
		uint32_t stage_minor = get<1>(*p);

		bool sensitive_mode = get<2>(*p);
		uint32_t prev_stage_major = get<0>(*(p - 1));
		uint32_t prev_stage_minor = get<1>(*(p - 1));

		uint32_t stage_id = (stage_major << 5) + stage_minor + (sensitive_mode ? (1 << 10) : 0);
#ifdef COLLECT_STATS
		objects.running_stats->Register(STAT_TIME_THR_BIN_READER_BASE + stage_id, "Thread time: bin reader " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_BIN_WRITER_BASE + stage_id, "Thread time: bin writer " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_MAPPING_CORE_BASE + stage_id, "Thread time: map. core  " + StageDesc(stage_id), running_stats_t::lists);
#endif

		objects.progress->SetComment("  (" + StageDesc(stage_id) + +")");

		watch.StartTimer();
		stage(stage_major, stage_minor, prev_stage_major, prev_stage_minor, mapping_mode_t::first, sensitive_mode, (p + 1) == stages_to_do.end());
		watch.StopTimer();
#ifdef COLLECT_STATS
		objects.running_stats->AddTotals(STAT_TIME_BASE + (stage_major << 5) + stage_minor + (sensitive_mode ? (1 << 10) : 0), watch.GetElapsedTime());
#endif
	}

	thr_rgw->join();
	delete objects.mp_map_res;
	delete objects.q_map_res;
	delete objects.mp_reads;
	delete thr_rgw;
	delete rgw;

	if (params.verbosity_level > 0)
	{
		cerr << "** End of mapping **      \n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
// Read bins and perform mapping (second stratum)
bool CMapper::reads_mapping_second_stratum()
{
	CStopWatch watch;

	string prefix_name = "main";

	if (params.verbosity_level > 0)
	{
		cerr << "***** Reads mapping *****\n";
		fflush(stderr);
	}

	uint32_t max_stage = params.max_no_errors;

#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_SORTING_TOTAL, "Reads sorting (total)", running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_RES_WRITER, "Thread time: results writer      ", running_stats_t::totals);
#endif

	// Prepare stage-substage descriptions
	vector<tuple<uint32_t, uint32_t, bool>> stages_to_do;
	for (uint32_t stage_major = 0; stage_major <= max_stage;)
	{
		for (uint32_t stage_minor = 0; stage_minor <= stage_major; ++stage_minor)
			stages_to_do.push_back(make_tuple(stage_major, stage_minor, false));

		stage_major = get_next_stage(stage_major, max_stage);
	}

	if (params.sensitive_mode)
		for (uint32_t stage_minor = 0; stage_minor <= max_stage; ++stage_minor)
			stages_to_do.push_back(make_tuple(max_stage, stage_minor, true));

	// Mapping results queue, thread, and memory pool
	uint32_t tot_substages = (uint32_t)stages_to_do.size() - 1;

	objects.progress->Init((int64_t) tot_substages * params.no_bins, true);

	objects.mp_reads = new CMemoryPool<uchar_t>(params.read_len * params.no_thr_mapping_cores, params.read_len, "mp_reads");
	objects.q_map_res = new CRegisteringQueue<res_group_t>(params.no_thr_mapping_cores * tot_substages);
	objects.mp_map_res = new CMemoryPool<uchar_t>(params.no_res_parts * params.res_group_size,
		params.res_group_size, "mp_map_res");
	CResultGroupsWriter *rgw = new CResultGroupsWriter(&params, &objects, prefix_name);
	thread *thr_rgw = new thread(ref(*rgw));

	for (auto p = stages_to_do.begin() + 1; p != stages_to_do.end(); ++p)
	{
		uint32_t stage_major = get<0>(*p);
		uint32_t stage_minor = get<1>(*p);
		bool sensitive_mode = get<2>(*p);
		uint32_t prev_stage_major = get<0>(*(p - 1));
		uint32_t prev_stage_minor = get<1>(*(p - 1));

		uint32_t stage_id = (stage_major << 5) + stage_minor + (sensitive_mode ? (1 << 10) : 0);
#ifdef COLLECT_STATS
		objects.running_stats->Register(STAT_TIME_THR_BIN_READER_BASE + stage_id, "Thread time: bin reader " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_BIN_WRITER_BASE + stage_id, "Thread time: bin writer " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_MAPPING_CORE_BASE + stage_id, "Thread time: map. core  " + StageDesc(stage_id), running_stats_t::lists);
#endif

		objects.progress->SetComment("  (" + StageDesc(stage_id) + +")");

		watch.StartTimer();
		stage(stage_major, stage_minor, prev_stage_major, prev_stage_minor, mapping_mode_t::second, sensitive_mode);
		watch.StopTimer();

#ifdef COLLECT_STATS
		objects.running_stats->AddTotals(STAT_TIME_BASE + (stage_major << 5) + stage_minor + (sensitive_mode ? (1 << 10) : 0), watch.GetElapsedTime());
#endif
	}

	thr_rgw->join();
	delete objects.mp_map_res;
	delete objects.q_map_res;
	delete objects.mp_reads;
	delete thr_rgw;
	delete rgw;

	if (params.verbosity_level > 0)
	{
		cerr << "** End of mapping **      \n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
// Read bins and perform mapping (all strata)
bool CMapper::reads_mapping_all_strata()
{
	CStopWatch watch;

	string prefix_name = "main";

	if (params.verbosity_level > 0)
	{
		cerr << "***** Reads mapping *****\n";
		fflush(stderr);
	}

	uint32_t max_stage = params.max_no_errors;

#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_SORTING_TOTAL, "Reads sorting (total)", running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_RES_WRITER, "Thread time: results writer      ", running_stats_t::totals);
#endif

	// Prepare stage-substage descriptions
	vector<pair<uint32_t, uint32_t>> stages_to_do;
	stages_to_do.push_back(make_pair(0, 0));
	for (uint32_t stage_minor = 0; stage_minor <= max_stage; ++stage_minor)
		stages_to_do.push_back(make_pair(max_stage, stage_minor));

	// Mapping results queue, thread, and memory pool
	uint32_t tot_substages = (uint32_t)stages_to_do.size() - 1;

	objects.mp_reads = new CMemoryPool<uchar_t>(params.read_len * params.no_thr_mapping_cores, params.read_len, "mp_reads");
	objects.q_map_res = new CRegisteringQueue<res_group_t>(params.no_thr_mapping_cores * tot_substages);
	objects.mp_map_res = new CMemoryPool<uchar_t>(params.no_res_parts * params.res_group_size,
		params.res_group_size, "mp_map_res");
	CResultGroupsWriter *rgw = new CResultGroupsWriter(&params, &objects, prefix_name);
	thread *thr_rgw = new thread(ref(*rgw));

	// Add here sensitive mode
	for (auto p = stages_to_do.begin() + 1; p != stages_to_do.end(); ++p)
	{
		uint32_t stage_major = p->first;
		uint32_t stage_minor = p->second;
		uint32_t prev_stage_major = (p - 1)->first;
		uint32_t prev_stage_minor = (p - 1)->second;

#ifdef COLLECT_STATS
		uint32_t stage_id = (stage_major << 5) + stage_minor;
		objects.running_stats->Register(STAT_TIME_THR_BIN_READER_BASE + stage_id, "Thread time: bin reader " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_BIN_WRITER_BASE + stage_id, "Thread time: bin writer " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_MAPPING_CORE_BASE + stage_id, "Thread time: map. core  " + StageDesc(stage_id), running_stats_t::lists);
#endif

		watch.StartTimer();
		stage(stage_major, stage_minor, prev_stage_major, prev_stage_minor, mapping_mode_t::all);
		watch.StopTimer();
		
#ifdef COLLECT_STATS
		objects.running_stats->AddTotals(STAT_TIME_BASE + (stage_major << 5) + stage_minor, watch.GetElapsedTime());
#endif
	}

	thr_rgw->join();
	delete objects.mp_map_res;
	delete objects.q_map_res;
	delete objects.mp_reads;
	delete thr_rgw;
	delete rgw;

	if (params.verbosity_level > 0)
	{
		cerr << "End of all strata mapping\n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
// Post processing of mapped reads
bool CMapper::reads_postprocessing(mapping_mode_t mapping_mode)
{
#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_TIME_THR_FASTQ_READER_PP, "Thread time: FASTQ reader post processing ", running_stats_t::lists);
	objects.running_stats->Register(STAT_TIME_THR_RES_READER, "Thread time: result groups readers        ", running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_SAM_GENERATOR, "Thread time: SAM generators               ", running_stats_t::lists);
	objects.running_stats->Register(STAT_TIME_THR_SAM_SORTING, "Thread time: SAM sorting                  ", running_stats_t::lists);
	objects.running_stats->Register(STAT_TIME_THR_SAM_PROCESSING, "Thread time: SAM processing               ", running_stats_t::lists);
#endif
	CStopWatch watch;
	watch.StartTimer();

	string prefix_name = "main";

	if (params.verbosity_level > 1)
	{
		cerr << "** Postprocessing (SAM generation) **\n";
		fflush(stderr);
	}

	// Prepare pool memory allocators
//	objects.mp_fastq_blocks = new CMemoryPool<uchar_t>(
//		(params.block_size + params.block_overhead_size) * (params.no_thr_fastq_reader + params.no_thr_sam_generators) * params.no_fastq_blocks_per_thread,
//		params.block_size + params.block_overhead_size);

	objects.mp_fastq_blocks = new CMemoryPool<uchar_t>((int64_t) 
		((params.block_size + params.block_overhead_size) * params.no_fastq_blocks_per_thread * (2.5 * params.no_thr_fastq_reader + params.no_thr_sam_generators)),
		params.block_size + params.block_overhead_size, "mp_fastq_blocks");
	objects.IncreaseMem("mp_fastq_blocks", (int64_t)((params.block_size + params.block_overhead_size) * params.no_fastq_blocks_per_thread *
		(2.5 * params.no_thr_fastq_reader + params.no_thr_sam_generators)));

	objects.mp_fastq_records = new CMemoryPool<uchar_t>((int64_t) params.max_fastq_rec_length * params.no_thr_sam_generators * 3 * 4, params.max_fastq_rec_length, 
		"mp_fastq_records");
	objects.IncreaseMem("mp_fastq_records", (int64_t) params.max_fastq_rec_length * params.no_thr_sam_generators * 3 * 4);

	if (params.gzipped_SAM_level == 0)
	{
		objects.mp_sam_parts = new CMemoryPool<uchar_t>((int64_t) params.sam_part_size * params.no_thr_sam_generators * 2, params.sam_part_size, "mp_sam_parts");
		objects.IncreaseMem("mp_sam_parts", (int64_t) params.sam_part_size * params.no_thr_sam_generators * 2);
	}
	else
	{
		objects.mp_sam_parts = new CMemoryPool<uchar_t>((int64_t) params.sam_part_size * params.no_thr_sam_generators * 3, params.sam_part_size, "mp_sam_parts");
		objects.IncreaseMem("mp_sam_parts", (int64_t) params.sam_part_size * params.no_thr_sam_generators * 3);
	}

	int mult_cigars = 4;

	objects.mp_ext_cigar = new CMemoryPool<uchar_t>((int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings, params.max_cigar_len, 
		"mp_ext_cigar");
	objects.IncreaseMem("mp_cigar_ext", (int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings);

	objects.mp_cigar = new CMemoryPool<uchar_t>((int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings, params.max_cigar_len, 
		"mp_cigar");
	objects.IncreaseMem("mp_cigar", (int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings);
	
	objects.mp_cigar_bin = new CMemoryPool<uint32_t>((int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings, params.max_cigar_len, 
		"mp_cigar_bin");
	objects.IncreaseMem("mp_cigar_bin", (int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings);

	objects.mp_mdz = new CMemoryPool<uchar_t>((int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings, params.max_cigar_len, 
		"mp_mdz");
	objects.IncreaseMem("mp_mdz", (int64_t) params.max_cigar_len * params.no_thr_sam_generators * mult_cigars * params.max_no_mappings);

	// Prepare queues
	objects.q_file_names = new CRegisteringQueue<file_name_no_t>(1);
	objects.q_res_ids = new CRegisteringQueue<uint32_t>(params.no_thr_fastq_reader);
	objects.q_map_reads = new CRegisteringQueue<mapped_reads_t>(1);
	objects.q_sam_blocks = new CRegisteringQueue<sam_block_t>(params.no_thr_sam_generators);

	objects.mem_monitor = new CMemoryMonitor(params.max_res_dev_memory);

	objects.ptr_pool = new CPtrPool((int) (params.no_thr_sam_generators * 2.5));

	// Load ID store
	CIDStore *id_store = new CIDStore(params.id_bits_total, params.id_bits_subgroup, params.id_bits_local, params.verbosity_level);
	shared_ptr<CMapperFile> mf(new CMapperFile(EXT_ID_STORE, MARKER_ID_STORE));
	mf->OpenRead(params.temp_prefix + PRJ_CFG_FILE);
	id_store->Load(mf);
	mf->Close();

	if (!params.keep_temporary_files)
		mf->Remove();

	params.no_res_groups = id_store->GetNumGroups();

	objects.progress->Init(params.no_res_groups, true);
	objects.progress->SetComment("");

	// Load input file names to the queue
	uint32_t file_no = 0;
	for (auto &p : params.input_file_names)
	{
		objects.q_file_names->Push(file_name_no_t(file_no, p.first, p.second));
		file_no += 2;
	}
	objects.q_file_names->MarkCompleted();

	// Construct joinger manager
	CJoinerMgr* joiner_mgr = new CJoinerMgr(&params, &objects);

	// Create results group deliverer thread
	CMappingResultsDeliverer *mrd = new CMappingResultsDeliverer(&params, &objects, prefix_name, joiner_mgr);
	thread *thr_mrd = new thread(ref(*mrd));

	load_ref_seq_desc();

	// Create SAM writer
	vector<string> header_SAM;
	vector<pair<string, uint32_t>> header_BAM;
	prepare_output_files(header_SAM, header_BAM);

	CSamWriter* sw = new CSamWriter(&params, &objects, header_SAM, header_BAM);
	thread *thr_sw = new thread(ref(*sw));

	// Create FASTQ reader threads
	vector<CPostFastqReader*> fqr;
	for (uint32_t i = 0; i < params.no_thr_fastq_reader; ++i)
		fqr.push_back(new CPostFastqReader(false, &params, &objects, id_store, joiner_mgr));
	vector<thread*> thr_fqr;
	for (uint32_t i = 0; i < params.no_thr_fastq_reader; ++i) {
		thr_fqr.push_back(new thread(ref(*fqr[i])));
	}

	// Create SAM generators
	vector<CSamGenerator*> sg;
	for (uint32_t i = 0; i < params.no_thr_sam_generators; ++i)
		sg.push_back(new CSamGenerator(&params, &objects, &ref_seq_desc));
	vector<thread*> thr_sg;
	for (uint32_t i = 0; i < params.no_thr_sam_generators; ++i) {
		thr_sg.push_back(new thread(ref(*sg[i])));
		//	(*sg[i])();
	}

	// Wait for completing the stage
	for (auto &p : thr_fqr)
		p->join();
	thr_mrd->join();

	for (auto &p : thr_sg)
		p->join();
	thr_sw->join();

	// Destruct objects
	delete id_store;
	delete joiner_mgr;
	delete thr_mrd;
	delete mrd;
	delete thr_sw;
	delete sw;

	for (uint32_t i = 0; i < params.no_thr_sam_generators; ++i)
	{
		delete thr_sg[i];
		delete sg[i];
	}

	for (uint32_t i = 0; i < params.no_thr_fastq_reader; ++i)
	{
		delete thr_fqr[i];
		delete fqr[i];
	}

	delete objects.mp_fastq_blocks;
	delete objects.mp_sam_parts;
	delete objects.mp_fastq_records;
	delete objects.mp_ext_cigar;
	delete objects.mp_cigar;
	delete objects.mp_cigar_bin;
	delete objects.mp_mdz;

	delete objects.q_file_names;
	delete objects.q_map_reads;
	delete objects.q_res_ids;
	delete objects.q_sam_blocks;

	delete objects.ptr_pool;

	delete objects.mem_monitor;

	cerr << "***** End of work *****\n";
	fflush(stderr);

	watch.StopTimer();
	//	objects.running_stats->AddTotals(STAT_TIME_SPLIT, watch.GetElapsedTime());

	return true;
}

// ************************************************************************************
// A single stage of mapping - only for development purposes!
// !!! Works correctly only for upto 4 errors
bool CMapper::reads_mapping_single_stage(uint32_t stage_major, uint32_t stage_minor, mapping_mode_t mapping_mode, bool sensitive_mode)
{
	CStopWatch watch;
	uint32_t prev_stage_major = 0;
	uint32_t prev_stage_minor = 0;

	string prefix_name = "main";

	if (params.verbosity_level > 0)
	{
		cerr << "***** Reads mapping - single stage *****\n";
		fflush(stderr);
	}

	//	uint32_t max_stage = max(params.max_no_mismatches, params.max_no_indels);
	//	max_stage = 1;

#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_SORTING_TOTAL, "Reads sorting (total)", running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_RES_WRITER, "Thread time: results writer      ", running_stats_t::totals);
#endif

	// Mapping results queue, thread, and memory pool
	objects.mp_reads = new CMemoryPool<uchar_t>(params.read_len * params.no_thr_mapping_cores, params.read_len, "mp_reads");
	objects.q_map_res = new CRegisteringQueue<res_group_t>(params.no_thr_mapping_cores);
	objects.mp_map_res = new CMemoryPool<uchar_t>(params.no_res_parts * params.res_group_size, params.res_group_size, "mp_map_res");
	CResultGroupsWriter *rgw = new CResultGroupsWriter(&params, &objects, prefix_name);
	thread *thr_rgw = new thread(ref(*rgw));

	if (!sensitive_mode)
	{
		if (stage_minor == 0)
		{
			prev_stage_major = stage_major - (stage_major == 4 ? 2 : 1);
			prev_stage_minor = prev_stage_major;
		}
		else
		{
			prev_stage_major = stage_major;
			prev_stage_minor = stage_minor - 1;
		}
	}
	else
	{
		if (stage_minor == 0)
		{
			prev_stage_major = stage_major - 1;
			prev_stage_minor = prev_stage_major;
		}
		else
		{
			prev_stage_major = stage_major;
			prev_stage_minor = stage_minor - 1;
		}
	}

#ifdef COLLECT_STATS
	uint32_t stage_id = (stage_major << 5) + stage_minor + (sensitive_mode ? (1 << 10) : 0);
	objects.running_stats->Register(STAT_TIME_THR_BIN_READER_BASE + stage_id, "Thread time: bin reader " + StageDesc(stage_id), running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_BIN_WRITER_BASE + stage_id, "Thread time: bin writer " + StageDesc(stage_id), running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_MAPPING_CORE_BASE + stage_id, "Thread time: map. core  " + StageDesc(stage_id), running_stats_t::lists);
#endif

	watch.StartTimer();
	stage(stage_major, stage_minor, prev_stage_major, prev_stage_minor, mapping_mode, sensitive_mode, false);
	watch.StopTimer();

#ifdef COLLECT_STATS
	objects.running_stats->AddTotals(STAT_TIME_BASE + (stage_major << 5) + stage_minor, watch.GetElapsedTime());
#endif

	prev_stage_minor = stage_minor;
	prev_stage_major = stage_major;

	thr_rgw->join();
	delete objects.mp_map_res;
	delete objects.q_map_res;
	delete objects.mp_reads;
	delete thr_rgw;
	delete rgw;

	if (params.verbosity_level > 0)
	{
		cerr << "End of single stage mapping\n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
//
bool CMapper::reads_mapping_stage_range(uint32_t stage_major_from, uint32_t stage_minor_from, uint32_t stage_major_to, uint32_t stage_minor_to,
	mapping_mode_t mapping_mode, bool sensitive_mode)
{
	CStopWatch watch;

	string prefix_name = "main";

	if (params.verbosity_level > 0)
	{
		cerr << "***** Reads mapping - range of stages *****\n";
		fflush(stderr);
	}

	uint32_t max_stage = params.max_no_errors;

#ifdef COLLECT_STATS
	objects.running_stats->Register(STAT_SORTING_TOTAL, "Reads sorting (total)", running_stats_t::totals);
	objects.running_stats->Register(STAT_TIME_THR_RES_WRITER, "Thread time: results writer      ", running_stats_t::totals);
#endif

	// Prepare stage-substage descriptions
	vector<pair<uint32_t, uint32_t>> stages_to_do;

	if (stage_minor_from > stage_major_from)
		stage_minor_from = stage_major_from;
	if (stage_minor_to > stage_major_to)
		stage_minor_to = stage_major_to;

	if (!sensitive_mode)
	{
		if (mapping_mode == mapping_mode_t::first)
		{
			if (stage_minor_from > 0)
				stages_to_do.push_back(make_pair(stage_major_from, stage_minor_from - 1));
			else if (stage_major_from > 0)
				stages_to_do.push_back(make_pair(stage_major_from - 1, stage_major_from - 1));
			else
				stages_to_do.push_back(make_pair(0, 0));
		}
		else if (mapping_mode == mapping_mode_t::all)
		{
			if (stage_minor_from > 0)
				stages_to_do.push_back(make_pair(stage_major_from, stage_minor_from - 1));
			else
				stages_to_do.push_back(make_pair(0, 0));
		}
	}
	else
	{
		if (stage_minor_from > 0)
			stages_to_do.push_back(make_pair(stage_major_from, stage_minor_from - 1));
		else
			stages_to_do.push_back(make_pair(stage_major_from, stage_major_from));
	}

	uint32_t stage_id_from = (stage_major_from << 5) + stage_minor_from;
	uint32_t stage_id_to = (stage_major_to << 5) + stage_minor_to;

	for (uint32_t stage_major = 0; stage_major <= max_stage;)
	{
		for (uint32_t stage_minor = 0; stage_minor <= stage_major; ++stage_minor)
		{
			uint32_t stage_id = (stage_major << 5) + stage_minor;
			if (stage_id >= stage_id_from && stage_id <= stage_id_to)
				stages_to_do.push_back(make_pair(stage_major, stage_minor));
		}
		/*		if(stage_major+2 == max_stage && max_stage >= 4)
		stage_major += 2;
		else*/
		stage_major += 1;
	}

	// Mapping results queue, thread, and memory pool
	uint32_t tot_substages = (uint32_t)stages_to_do.size() - 1;

	objects.mp_reads = new CMemoryPool<uchar_t>(params.read_len * params.no_thr_mapping_cores, params.read_len, "mp_reads");
	objects.q_map_res = new CRegisteringQueue<res_group_t>(params.no_thr_mapping_cores * tot_substages);
	objects.mp_map_res = new CMemoryPool<uchar_t>(params.no_res_parts * params.res_group_size,
		params.res_group_size, "mp_map_res");
	CResultGroupsWriter *rgw = new CResultGroupsWriter(&params, &objects, prefix_name);
	thread *thr_rgw = new thread(ref(*rgw));

	for (auto p = stages_to_do.begin() + 1; p != stages_to_do.end(); ++p)
	{
		uint32_t stage_major = p->first;
		uint32_t stage_minor = p->second;
		uint32_t prev_stage_major = (p - 1)->first;
		uint32_t prev_stage_minor = (p - 1)->second;

#ifdef COLLECT_STATS
		uint32_t stage_id = (stage_major << 5) + stage_minor;
		objects.running_stats->Register(STAT_TIME_THR_BIN_READER_BASE + stage_id, "Thread time: bin reader " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_BIN_WRITER_BASE + stage_id, "Thread time: bin writer " + StageDesc(stage_id), running_stats_t::totals);
		objects.running_stats->Register(STAT_TIME_THR_MAPPING_CORE_BASE + stage_id, "Thread time: map. core  " + StageDesc(stage_id), running_stats_t::lists);
#endif

		watch.StartTimer();
		stage(stage_major, stage_minor, prev_stage_major, prev_stage_minor, mapping_mode, sensitive_mode);
		watch.StopTimer();
		
#ifdef COLLECT_STATS
		objects.running_stats->AddTotals(STAT_TIME_BASE + (stage_major << 5) + stage_minor, watch.GetElapsedTime());
#endif
	}

	thr_rgw->join();
	delete objects.mp_map_res;
	delete objects.q_map_res;
	delete objects.mp_reads;
	delete thr_rgw;
	delete rgw;

	if (params.verbosity_level > 1)
	{
		cerr << "End of stage range mapping\n";
		fflush(stderr);
	}

	return true;
}

// ************************************************************************************
bool CMapper::stage(uint32_t stage_major, uint32_t stage_minor, uint32_t prev_stage_major, uint32_t prev_stage_minor, mapping_mode_t mapping_mode, bool sensitive_mode, bool final_substage)
{
	string prefix_name_prev;
	string prefix_name_curr;
	uint32_t stage_id = (stage_major << 5) + stage_minor + (sensitive_mode ? (1 << 10) : 0);

	prefix_name_prev = "stage_" + Int2String(prev_stage_major) + (sensitive_mode && stage_minor > 0 ? "s_" : "_") + Int2String(prev_stage_minor) + "_";
	prefix_name_curr = "stage_" + Int2String(stage_major) + (sensitive_mode ? "s_" : "_") + Int2String(stage_minor) + "_";

	if (params.verbosity_level > 1)
	{
		cerr << "--------------------------\n";
		cerr << "Stage: " << stage_major << ":" << stage_minor << (sensitive_mode ? " sensitive\n" : "\n");
		fflush(stderr);
	}

	objects.q_bins_read = new CRegisteringQueue<reads_bin_t>(1);
	objects.q_bins_write = new CRegisteringQueue<reads_bin_t>(params.no_thr_mapping_cores);

	objects.mem_monitor = new CMemoryMonitor(params.max_mapping_memory);
	objects.mp_bins_write = new CMemoryPool<uchar_t>(params.no_parts * params.bin_size, params.bin_size, "mp_bins_write");

	// Create bin writer thread
	CBinsWriter *bw = new CBinsWriter(&params, &objects, prefix_name_curr, stage_id, final_substage);
	thread *thr_bw = new thread(ref(*bw));

	// Mapping core threads
	vector<CMappingCore*> mc(params.no_thr_mapping_cores);
	vector<thread*> thr_mc(params.no_thr_mapping_cores);
	for (uint32_t i = 0; i < params.no_thr_mapping_cores; ++i)
	{
		mc[i] = new CMappingCore(&params, &objects, stage_id, mapping_mode);
		thr_mc[i] = new thread(ref(*mc[i]));
	}

	CBinsReader *br = new CBinsReader(&params, &objects, prefix_name_prev, stage_id);
	thread *thr_br = new thread(ref(*br));

	thr_br->join();
	for (auto &p : thr_mc)
		p->join();
	thr_bw->join();

	delete thr_br;
	delete br;
	delete thr_bw;
	delete bw;

	for (uint32_t i = 0; i < params.no_thr_mapping_cores; ++i)
	{
		delete thr_mc[i];
		delete mc[i];
	}

	delete objects.q_bins_read;
	delete objects.q_bins_write;
	delete objects.mem_monitor;
	delete objects.mp_bins_write;

	if (params.verbosity_level > 1)
		cerr << "\n";

	return true;
}

// ************************************************************************************
// Prepare SAM files
bool CMapper::prepare_output_files(vector<string> &header_SAM, vector<pair<string, uint32_t>> &header_BAM)
{
	header_SAM.clear();
	header_BAM.clear();

	// Mapped reads
	header_SAM.push_back("@HD\tVN:1.3\tSO:unsorted\n");

	vector<seq_desc_t> seq_desc;
	ref_seq_desc.GetDescription(seq_desc);

	for (auto &p : seq_desc)
	{
		header_SAM.push_back("@SQ\t");
		header_SAM.push_back("SN:" + p.name + "\t");
		header_SAM.push_back("LN:" + Int2String(p.size + p.no_initial_Ns + p.no_final_Ns));
		header_SAM.push_back("\n");

		header_BAM.push_back(make_pair(p.name, (uint32_t) (p.size + p.no_initial_Ns + p.no_final_Ns)));
	}

	string cmd = params.command_line;
	for (auto c = cmd.begin(); c != cmd.end(); ++c)
		if (*c == '\"')
			*c = ' ';
		else if (*c == '@')
			*c = ' ';

	if (params.read_group_line.length() > 0) {
		header_SAM.push_back(params.read_group_line + "\n");
	}

	header_SAM.push_back("@PG\tID:" + string(MAPPER_ID));
	header_SAM.push_back("\tPN:" + string(MAPPER_NAME));
	header_SAM.push_back("\tVN:" + string(MAPPER_VERSION));
	header_SAM.push_back("\tCL:" + cmd);
	header_SAM.push_back("\n");
	

	return true;
}

// ************************************************************************************
bool CMapper::load_ref_seq_desc()
{
	return ref_seq_desc.Load(params.index_name);
}

// ************************************************************************************
bool CMapper::check_input_files()
{
	total_input_file_size = 0;

	for (auto &x : params.input_file_names)
	{
		for (auto &p : { x.first, x.second })
		{
			int64_t size = 0;
			if (p != "")
				size = FileSize(p);
			if (size < 0)
			{
				cerr << "No file: " << x.first << endl;
				return false;
			}

			total_input_file_size += size;
		}
	}

	return true;
}

// EOF 