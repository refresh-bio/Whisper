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


#include <cstdio>
#include "mapper.h"
#include <algorithm>

#include <iostream>
#include <fstream>

using namespace std;

uint32_t stage_major = 0xffff;
uint32_t stage_minor = 0;
uint32_t stage_major_from = 0xffff;
uint32_t stage_minor_from = 0;
uint32_t stage_major_to   = 0xffff;
uint32_t stage_minor_to   = 0;

CCmdParams cmd_params;

void usage();
bool parse_parameters(int argc, char **argv);
bool run_mapper(int argc, char **argv);

// ************************************************************************************
void usage()
{
	cerr << "Whisper v. " << MAPPER_VERSION << "\n";
	cerr << "Usage:\n";
	cerr << "   whisper [options] <index_name> @<files> \n";
	cerr << "   whisper [options] <index_name> file_se \n";
	cerr << "   whisper [options] <index_name> file_pe1 file_pe2\n";
	cerr << "Parameters:\n";
	cerr << "  index_name   - name of the index (as created by asm_pp)\n";
	cerr << "  files        - name of the file containing list of FASTQ files with seq. reads\n";
	cerr << "  file_se      - FASTQ file (single-end)\n";
	cerr << "  file_pe[1|2] - FASTQ files (paired-end)\n";
	cerr << "Options:\n";
	cerr << "  -b <value> - no. of temporary files (minimum: 100, default: " << cmd_params.no_bins << ")\n";
	cerr << "  -d[fr/ff/rf] - mapping orientation (default: " <<
		(cmd_params.mapping_orientation == mapping_orientation_t::forward_forward ? "-dff" :
			cmd_params.mapping_orientation == mapping_orientation_t::forward_reverse ? "-dfr (forward - reverse)" : 
			"-drf") << "\n";
#ifdef _DEVELOPMENT_MODE
	cerr << "  -dev-mode - turn on developer mode (default: " << cmd_params.developer_mode << ")\n";
#endif
	cerr << "  -dist_paired <value> - max. distance for paired read (default: " << cmd_params.max_mate_distance << ")\n";
	cerr << "  -e <value> - max. no of errors (default: " << cmd_params.max_no_errors << ", max: " << 100.0 * largest_frac_errors << "% of read length)\n";
	cerr << "  -e-paired <value> - max. fraction of errors in paired read (default: " << cmd_params.max_mate_edit_distance_frac << ")\n";
	cerr << "  -enable-boundary-clipping <value> - enable clipping at boundaries when a lot of mismatches appears (default: " << cmd_params.enable_boundary_clipping << ")\n";
	cerr << "  -filter <value> - store only mappings for given chromosome (default: " << cmd_params.filter << ")\n";
	cerr << "  -gap-open <value> - score for gap open (default: " << cmd_params.gap_open << ")\n";
	cerr << "  -gap-extend <value> - score for gap extend (default: " << cmd_params.gap_extend << ")\n";
	cerr << "  -gzipped-SAM-level <value> - gzip compression level of SAM/BAM, 0 - no compression (default: " << cmd_params.gzipped_SAM_level << ")\n";
	cerr << "  -hit-merging-threshold <value> - minimal distance between different mappings (default: " << cmd_params.hit_merging_threshold << ")\n";
	cerr << "  -high-confidence-sigmas <value> - (default: " << cmd_params.high_confidence_sigmas << ")\n";
	cerr << "  -hit-merging-wrt-first <value> - calculate distance in marged group w.r.t. first (default: " << cmd_params.hit_merging_wrt_first << ")\n";
	cerr << "  -m[f/s/a] - mode: first stratum/second stratum/all strata (default: " <<
//	cerr << "  -m[f|a] - mode: first stratum/all strata (default: " <<
		(cmd_params.mapping_mode == mapping_mode_t::first ? "first stratum" :
			cmd_params.mapping_mode == mapping_mode_t::second ? "second stratum" : 
			"all strata") << ")\n";
	cerr << "  -mask-lqb <value> - mask bases of quality lower than value (default: " << cmd_params.mask_low_quality_bases << ")\n";
//	cerr << "  -o <value> - approximate amount of RAM in GB (default: " << cmd_params.max_total_memory << "; must be at least equal no. threads, but not less than 8)\n";
	cerr << "  -out <name> - name of the output file (default: " << cmd_params.project_name << ")\n";
	cerr << "  -penalty-saturation <value> - no. of sigmas for max. penalty in matching pairs (default: " << cmd_params.penalty_saturation_sigmas << ")\n";
	cerr << "  -rg <read_group> - complete read group header line, ’\t’ character will be converted to a TAB in the output SAM while the read group ID will be attached to every read (example line: ’@RG\tID:foo\tSM:bar’)\n";
	cerr << "  -r[s|p] - single or paired-end reads (default: " << (cmd_params.paired_reads ? "paired" : "single") << ")\n";
	cerr << "  -score-discretization-threshold (default: " << cmd_params.score_discretization_threshold << ")\n";
	cerr << "  -score-match <value> - score for matching symbol (default: " << cmd_params.match_score << ")\n";
	cerr << "  -score-clipping <value> score for clipping (default: " << cmd_params.clipping_score << ")\n";
	cerr << "  -score-mismatch <value> - score for mismatching symbol (default: " << cmd_params.mismatch_score << ")\n";
	cerr << "  -sens <value> - turn on/off sensitive mode (default: " << (bool) cmd_params.sensitive_mode << ")\n";
	cerr << "  -sens-factor <value> - sensitivity factor (default: " << cmd_params.sensitivity_factor << ")\n";
	cerr << "  -stdout - use stdout to store the output (default: " << cmd_params.use_stdout << ")\n";
	cerr << "  -store-BAM - turn on saving in BAM (default: " << (bool)cmd_params.store_BAM << ")\n";
	cerr << "  -t <value> - no. of threads (0-adjust to hardware) (default: " << cmd_params.no_threads << ")\n";
	cerr << "  -temp <name> - prefix for temporary files (default: " << cmd_params.temp_prefix << ")\n";
	cerr << "  -x - load complete suffix arrays in main memory (default: " << cmd_params.sa_in_ram << ")\n";

#ifdef _DEVELOPMENT_MODE
	cerr << "Developer-mode options:\n";
	cerr << "  -log-stats - turn on logging statistics (default: " << cmd_params.log_stats << ")\n";
	cerr << "  -keep <value> - keep temporary files (default: " << cmd_params.keep_temporary_files << ")\n";
	cerr << "  -mapq-div <value> - divisor in MAPQ calculation formula (default: " << cmd_params.mapq_div << ")\n";
	cerr << "  -mapq-mult <value> - multiplier in MAPQ calculation formula (default: " << cmd_params.mapq_mult << ")\n";
	cerr << "  -mrl <value> - minimal read length (default: " << cmd_params.min_read_len << ")\n";
	cerr << "  -Mrl <value> - maximal read length (default: " << cmd_params.max_read_len << ")\n";
	cerr << "  -s <major>:<minor> - perform only stage <major>:<minor>\n";
	cerr << "  -s x - perform only postprocessing stage\n";
	cerr << "  -sr <major1>:<minor1>-<major2>:<minor2> - perform only stages from <major1>:<minor1> to <major2>:<minor2>\n";
	cerr << "  -v <value> - verbosity level (default: " << cmd_params.verbosity_level << ")\n";
#endif
	cerr << "Examples:\n";
	cerr << "  whisper human @files\n";
	cerr << "  whisper -temp temp/ human reads1.fq reads2.fq\n";
	cerr << "  whisper -out result.sam -temp temp/ -t 12 human reads1.fq reads2.fq\n";
}

// ************************************************************************************
// Parse the parameters
bool parse_parameters(int argc, char **argv)
{
	int i;
#ifdef _DEVELOPMENT_MODE
	uint64_t tmp;
#endif

	if(argc < 6)
		return false;

#ifndef _DEVELOPMENT_MODE
	cmd_params.verbosity_level = 1;
#endif

	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;
		// Number of threads
		if(strcmp(argv[i], "-t") == 0 && i + 1 < argc)
			cmd_params.no_threads = atoi(argv[++i]);
		// Approximate memory limit
/*		else if(strcmp(argv[i], "-o") == 0 && i + 1 < argc)
		{
			double t = atof(argv[++i]);
			tmp = ((uint64_t) t) << 30;
			cmd_params.max_total_memory = max(tmp, (uint64_t) (8ull << 30));
		}*/
		// Mapping orientation - forward-reverse
		else if(strcmp(argv[i], "-dfr") == 0)
			cmd_params.mapping_orientation = mapping_orientation_t::forward_reverse;
		// Mapping orientation - forward-forward
		else if(strcmp(argv[i], "-dff") == 0)
			cmd_params.mapping_orientation = mapping_orientation_t::forward_forward;
		// Mapping orientation - reverse-forward
		else if(strcmp(argv[i], "-drf") == 0)
			cmd_params.mapping_orientation = mapping_orientation_t::reverse_forward;
		// FASTA input files
		else if(strcmp(argv[i], "-fa") == 0)
			cmd_params.is_fasta = true;
		// FASTQ input files
		else if(strcmp(argv[i], "-fq") == 0)
			cmd_params.is_fasta = false;
		else if (strcmp(argv[i], "-dev-mode") == 0)
			cmd_params.developer_mode = true;
		else if (strcmp(argv[i], "-log-stats") == 0)
			cmd_params.log_stats = true;
		else if (strcmp(argv[i], "-stdout") == 0)
			cmd_params.use_stdout = true;
		// Gzipped SAM
		else if (strcmp(argv[i], "-gzipped-SAM") == 0 && i + 1 < argc)
			cmd_params.gzipped_SAM_level = NormalizeValue(atoi(argv[++i]), 0, 9);
		else if (strcmp(argv[i], "-filter") == 0 && i + 1 < argc)
			cmd_params.filter = string(argv[++i]);
		else if (strcmp(argv[i], "-e") == 0 && i + 1 < argc)
			cmd_params.max_no_errors = atoi(argv[++i]);
		else if (strcmp(argv[i], "-e-paired") == 0 && i + 1 < argc)
			cmd_params.max_mate_edit_distance_frac = atof(argv[++i]);
		else if (strcmp(argv[i], "-dist-paired") == 0 && i + 1 < argc)
			cmd_params.max_mate_distance = atoi(argv[++i]);
		else if (strcmp(argv[i], "-sens") == 0 && i + 1 < argc)
			cmd_params.sensitive_mode = atoi(argv[++i]);
		else if (strcmp(argv[i], "-sens-factor") == 0 && i + 1 < argc)
			cmd_params.sensitivity_factor = atof(argv[++i]);
		else if (strcmp(argv[i], "-store-BAM") == 0)
			cmd_params.store_BAM = true;
		else if (strcmp(argv[i], "-out") == 0 && i + 1 < argc)
			cmd_params.project_name = string(argv[++i]);
		else if (strcmp(argv[i], "-temp") == 0 && i + 1 < argc)
			cmd_params.temp_prefix = string(argv[++i]);
		else if (strcmp(argv[i], "-rg") == 0 && i + 1 < argc) 
			cmd_params.read_group_line = string(argv[++i]);
#ifdef _DEVELOPMENT_MODE
		else if(strcmp(argv[i], "-v") == 0 && i + 1 < argc)
			cmd_params.verbosity_level = NormalizeValue(atoi(argv[++i]), 1, 4);
		else if (strcmp(argv[i], "-keep") == 0 && i + 1 < argc)
			cmd_params.keep_temporary_files = atoi(argv[++i]);
		else if(strcmp(argv[i], "-mrl") == 0 && i + 1 < argc)
			cmd_params.min_read_len = atoi(argv[++i]);
		else if(strcmp(argv[i], "-Mrl") == 0 && i + 1 < argc)
			cmd_params.max_read_len = atoi(argv[++i]);
		else if(strcmp(argv[i], "-sr") == 0 && i + 1 < argc)
		{
			if(strlen(argv[++i]) == 7)
			{
				stage_major_from = argv[i][0] - '0';
				stage_minor_from = argv[i][2] - '0';
				stage_major_to   = argv[i][4] - '0';
				stage_minor_to   = argv[i][6] - '0';
			}
		}
		else if(strcmp(argv[i], "-s") == 0 && i + 1 < argc)
		{
			++i;
			if(strlen(argv[i]) == 3)
			{
				stage_major = argv[i][0] - '0';
				stage_minor = argv[i][2] - '0';
			}
			if(strlen(argv[i]) == 1 && argv[i][0] == 'x')
				stage_major = stage_minor = 0xff;
		}
#endif
		// Mapping mode
		else if (strcmp(argv[i], "-mf") == 0)
			cmd_params.mapping_mode = mapping_mode_t::first;
		else if (strcmp(argv[i], "-ms") == 0)
			cmd_params.mapping_mode = mapping_mode_t::second;
		else if (strcmp(argv[i], "-ma") == 0)
			cmd_params.mapping_mode = mapping_mode_t::all;
		else if (strcmp(argv[i], "-mask-lqb") == 0 && i + 1 < argc)
			cmd_params.mask_low_quality_bases = atoi(argv[++i]);
		else if (strcmp(argv[i], "-penalty-saturation") == 0 && i + 1 < argc)
			cmd_params.penalty_saturation_sigmas = atof(argv[++i]);
		// Single or paired-end reads
		else if(strcmp(argv[i], "-rs") == 0)
			cmd_params.paired_reads = false;
		else if (strcmp(argv[i], "-rp") == 0)
			cmd_params.paired_reads = true;
		else if (strcmp(argv[i], "-hit-merging-threshold") == 0 && i + 1 < argc)
			cmd_params.hit_merging_threshold = atoi(argv[++i]);
		else if (strcmp(argv[i], "-hit-merging-wrt-first") == 0 && i + 1 < argc)
			cmd_params.hit_merging_wrt_first = (bool)atoi(argv[++i]);
		else if (strcmp(argv[i], "-high-confidence-sigmas") == 0 && i + 1 < argc)
			cmd_params.high_confidence_sigmas = atoi(argv[++i]);		
		else if (strcmp(argv[i], "-score-match") == 0 && i + 1 < argc)
			cmd_params.match_score = atof(argv[++i]);
		else if (strcmp(argv[i], "-score-mismatch") == 0 && i + 1 < argc)
			cmd_params.mismatch_score = atof(argv[++i]);
		else if (strcmp(argv[i], "-score-clipping") == 0 && i + 1 < argc)
			cmd_params.clipping_score = atof(argv[++i]);
		else if (strcmp(argv[i], "-score-discretization-threshold") == 0 && i + 1 < argc)
			cmd_params.score_discretization_threshold = atof(argv[++i]);
		else if (strcmp(argv[i], "-gap-open") == 0 && i + 1 < argc)
			cmd_params.gap_open = atof(argv[++i]);
		else if (strcmp(argv[i], "-gap-extend") == 0 && i + 1 < argc)
			cmd_params.gap_extend = atof(argv[++i]);
		else if (strcmp(argv[i], "-enable-boundary-clipping") == 0 && i + 1 < argc)
			cmd_params.enable_boundary_clipping = (bool) atoi(argv[++i]);
		else if (strcmp(argv[i], "-mapq-mult") == 0 && i + 1 < argc)
			cmd_params.mapq_mult = atof(argv[++i]);
		else if (strcmp(argv[i], "-mapq-div") == 0 && i + 1 < argc)
			cmd_params.mapq_div = atof(argv[++i]);

		// No. of temporary bins
		else if(strcmp(argv[i], "-b") == 0)
			cmd_params.no_bins = NormalizeValue(atoi(argv[++i]), 100, 2040);
		else if(strcmp(argv[i], "-x") == 0)
			cmd_params.sa_in_ram = true;
	}

	if(argc - i < 2)
		return false;

	cmd_params.index_name = string(argv[i++]);

	cmd_params.input_file_names.clear();

	// Single file - SE
	if (i + 1 == argc && argv[i][0] != '@')
		cmd_params.input_file_names.push_back(make_pair(string(argv[i]), ""));
	// Pair of files - PE
	else if(i + 2 == argc)
		cmd_params.input_file_names.push_back(make_pair(string(argv[i]), string(argv[i+1])));
	// File with list of files
	else
	{
		string input_file_name = string(argv[i++]);
	
		ifstream in(input_file_name.c_str()+1);
		if(!in.good())
		{
			cerr << "Error: No " << input_file_name.c_str()+1 << " file\n";
			return false;
		}

		string s, fn1, fn2;
		while(getline(in, s))
			if(s != "")
			{
				size_t pos = s.find(';');
				if(pos == string::npos)
					cmd_params.input_file_names.push_back(make_pair(s, ""));
				else
				{
					fn1 = s.substr(0, pos);
					fn2 = s.substr(pos+1, string::npos);
					cmd_params.input_file_names.push_back(make_pair(fn1, fn2));
				}
			}

		in.close();
		random_shuffle(cmd_params.input_file_names.begin(), cmd_params.input_file_names.end());
	}

	// Prepare command line for SAM
	cmd_params.command_line = "\"";

	for (int i = 0; i < argc; ++i)
	{
		if (i > 0)
			cmd_params.command_line.append(" ");
		cmd_params.command_line.append(argv[i]);
	}
	cmd_params.command_line.append("\"");
	
	// Set memory size
	if (cmd_params.no_threads < 8)
		cmd_params.max_total_memory = 8ull << 30;
	else
		cmd_params.max_total_memory = ((uint64_t)cmd_params.no_threads) << 30;
	
	// Check memory settings
/*	if ((cmd_params.max_total_memory >> 30) < cmd_params.no_threads)
	{
		cerr << "Warning: Max. memory usage (in GB) must be at least equal to no. of threads (but not less than 8GB).\nMemory setting is corrected to " << cmd_params.no_threads << "GB\n";
		cmd_params.max_total_memory = ((uint64_t) cmd_params.no_threads) << 30;
	}*/

	return true;
}

// ************************************************************************************
bool run_mapper(int argc, char **argv)
{
	CStopWatch timer;
	timer.StartTimer();

	shared_ptr<CMapper> mapper(new CMapper);

	mapper->SetParams(cmd_params);

	CRunningStats *rs = new CRunningStats();

	mapper->SetRunningStats(rs);

	bool res;
	
	if(stage_major == 0xff)
		res = mapper->StartPostProcessing();
	else if(stage_major != 0xffff)
		res = mapper->StartMapping(stage_major, stage_minor, cmd_params.sensitive_mode);
	else if(stage_major_from != 0xffff)
		res = mapper->StartMapping(stage_major_from, stage_minor_from, stage_major_to, stage_minor_to, cmd_params.sensitive_mode);
	else
		res = mapper->StartMapping();

	timer.StopTimer();

	cerr << "Computation time: " << timer.GetElapsedTime() << " s\n";

	if (cmd_params.log_stats)
	{
		ofstream stream("whisper.stats_log");
		rs->CreateReport(stream);
	}

	delete rs;

	return res;
}

// ************************************************************************************
int main(int argc, char **argv)
{
#ifdef WIN32
	_setmaxstdio(2040);
#endif

	if(!parse_parameters(argc, argv))
	{
		usage();
		return 0;
	}

	run_mapper(argc, argv);

	return 0;
}

// EOF
