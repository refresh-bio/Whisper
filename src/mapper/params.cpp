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


#include "params.h"

ostream& operator<<(ostream &out, const CParams &p)
{
	out << "*** Params: \n";
	out << "bin_size                      : " << p.bin_size << endl;
	out << "block_overhead_size           : " << p.block_overhead_size << endl;
	out << "block_size                    : " << p.block_size << endl;
	out << "clipping_score                : " << p.clipping_score << endl;
	out << "command_line                  : " << p.command_line << endl;
	out << "constant_read_len             : " << p.constant_read_len << endl;
	out << "developer_mode                : " << p.developer_mode << endl;
	out << "enable_boundary_clipping      : " << p.enable_boundary_clipping << endl;
	out << "enable_paired_clipping		  : " << p.enable_paired_clipping << endl;
//	out << "error_score                   : " << p.error_score << endl;
	out << "filter                        : " << p.filter << endl;
//	out << "gap_open                      : " << p.gap_open << endl;
//	out << "gap_extend                    : " << p.gap_extend << endl;
	out << "gap_ins_open                  : " << p.gap_ins_open << endl;
	out << "gap_ins_extend                : " << p.gap_ins_extend << endl;
	out << "gap_del_open                  : " << p.gap_del_open << endl;
	out << "gap_del_extend                : " << p.gap_del_extend << endl;
	out << "gzipped_SAM_level             : " << p.gzipped_SAM_level << endl;
	out << "high_confidence_sigmas        : " << p.high_confidence_sigmas << endl;
	out << "hit_merging_threshold         : " << p.hit_merging_threshold << endl;
	out << "hit_merging_wrt_first         : " << p.hit_merging_wrt_first << endl;
	out << "id_bits_local                 : " << p.id_bits_local << endl;
	out << "id_bits_subgroup              : " << p.id_bits_subgroup << endl;
	out << "id_bits_total                 : " << p.id_bits_total << endl;
	out << "index_name                    : " << p.index_name << endl;

	out << "input_file_names              : [";
	for (auto x : p.input_file_names)
		out << x.first << "," << x.second << " ; ";
	out << "]" << endl;

	out << "io_calls_in_serial_mode       : " << p.io_calls_in_serial_mode << endl;
	out << "is_fasta                      : " << p.is_fasta << endl;
	out << "keep_temporary_files          : " << p.keep_temporary_files << endl;
	out << "log_stats                     : " << p.log_stats << endl;
	out << "mapping_counter_size          : " << p.mapping_counter_size << endl;
	out << "mapping_mode                  : " << (int)p.mapping_mode << endl;
	out << "mapping_orientation           : " << (int)p.mapping_orientation << endl;
	out << "mapq_div                      : " << p.mapq_div << endl;
	out << "mapq_mult                     : " << p.mapq_mult << endl;
	out << "mask_low_quality_bases        : " << p.mask_low_quality_bases << endl;
	out << "match_score                   : " << p.match_score << endl;
	out << "max_approx_indel_len          : " << p.max_approx_indel_len << endl;
	out << "max_bin_part_size             : " << p.max_bin_part_size << endl;
	out << "max_cigar_len                 : " << p.max_cigar_len << endl;
	out << "max_fastq_rec_length          : " << p.max_fastq_rec_length << endl;
	out << "max_no_errors                 : " << p.max_no_errors << endl;
	out << "max_group_part_size           : " << p.max_group_part_size << endl;
	out << "max_mapping_memory            : " << p.max_mapping_memory << endl;
	out << "max_mate_distance             : " << p.max_mate_distance << endl;
	out << "max_mate_edit_distance_frac   : " << p.max_mate_edit_distance_frac << endl;
	out << "max_no_errors                 : " << p.max_no_errors << endl;
	out << "max_no_mappings               : " << p.max_no_mappings << endl;
	out << "max_read_len                  : " << p.max_read_len << endl;
	out << "max_reads_compression         : " << p.max_reads_compression << endl;
	out << "max_res_dev_memory            : " << p.max_res_dev_memory << endl;
	out << "max_total_memory              : " << p.max_total_memory << endl;
	out << "min_approx_indel_len          : " << p.min_approx_indel_len << endl;
	out << "min_no_free_group_parts       : " << p.min_no_free_group_parts << endl;
	out << "min_no_free_parts             : " << p.min_no_free_parts << endl;
	out << "min_read_len                  : " << p.min_read_len << endl;
	out << "mismatch_score                : " << p.mismatch_score << endl;
	out << "no_bins                       : " << p.no_bins << endl;
	out << "no_fastq_blocks_per_thread    : " << p.no_fastq_blocks_per_thread << endl;
	out << "no_parts                      : " << p.no_parts << endl;
	out << "no_res_groups                 : " << p.no_res_groups << endl;
	out << "no_res_parts                  : " << p.no_res_parts << endl;
	out << "no_thr_fastq_reader           : " << p.no_thr_fastq_reader << endl;
	out << "no_thr_mapping_cores          : " << p.no_thr_mapping_cores << endl;
	out << "no_thr_sam_generators         : " << p.no_thr_sam_generators << endl;
	out << "no_thr_splitter               : " << p.no_thr_splitter << endl;
	out << "no_threads                    : " << p.no_threads << endl;
	out << "paired_reads                  : " << p.paired_reads << endl;
	out << "penalty_saturation_sigmas     : " << p.penalty_saturation_sigmas << endl;
	out << "project_name                  : " << p.project_name << endl;
	out << "read_len                      : " << p.read_len << endl;
	out << "res_group_size                : " << p.res_group_size << endl;
	out << "read_group_line               : " << p.read_group_line << endl;
	out << "sa_in_ram                     : " << p.sa_in_ram << endl;
	out << "sa_prefix_overhead            : " << p.sa_prefix_overhead << endl;
	out << "sam_buffer_memory             : " << p.sam_buffer_memory << endl;
	out << "sam_part_size                 : " << p.sam_part_size << endl;
	out << "score_discretization_threshold: " << p.score_discretization_threshold << endl;
	out << "sensitive_mode                : " << p.sensitive_mode << endl;
	out << "sensitivity_factor            : " << p.sensitivity_factor << endl;
	out << "ssd_mode                      : " << p.ssd_mode << endl;
	out << "temp_prefix                   : " << p.temp_prefix << endl;
	out << "verbosity_level               : " << p.verbosity_level << endl;

	out << endl;

	return out;
}