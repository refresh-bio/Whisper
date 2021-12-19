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


#ifndef _DEFS_H
#define _DEFS_H

#ifndef WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#endif

//#define _DEVELOPMENT_MODE

//#define COLLECT_STATS

// Old code for mismatch-only processing 
//#define MISMATCHES_CODE
//#define CONST_LEN_CODE

//#define ENABLE_VCF_VARIANTS

#include <string>
using namespace std;

#define MAPPER_VERSION		"2.0.2 (2021-12-19)"
#define MAPPER_ID			"Whisper"
#define MAPPER_NAME			"Whisper"


#define BOUNDARY_N_LEN		10240
#define SEPARATE_N_LEN		1536

#ifdef WIN32
typedef unsigned long long uint64_t;
typedef long long int64_t;
typedef unsigned int uint32_t;
typedef int int32_t;
#endif
typedef unsigned char uchar_t;

typedef uint32_t ref_pos_t;

const char sym_code_A      = 0;
const char sym_code_C      = 1;
const char sym_code_G      = 2;
const char sym_code_T	   = 3;
const char sym_code_N_read = 4;
const char sym_code_N_ref  = 5;

const uint8_t mask_code_A = (uint8_t)(1u << sym_code_A);
const uint8_t mask_code_C = (uint8_t)(1u << sym_code_C);
const uint8_t mask_code_G = (uint8_t)(1u << sym_code_G);
const uint8_t mask_code_T = (uint8_t)(1u << sym_code_T);
const uint8_t mask_code_N_read = 0;
const uint8_t mask_code_N_ref = 0;
const uint8_t mask_code_N = 0;

const uint32_t lut_short_prefix_len = 11;
const uint32_t lut_long_prefix_len  = 12;

const uint32_t bin_total_bits     = 40;
const uint32_t bin_sub_block_bits = 3;
const uint32_t bin_in_block_bits  = 21;

const uint64_t fastq_block_size = (1 << 25);

enum class genome_t {direct, rev_comp};
enum class mapping_mode_t {first, second, all};
enum class mapping_orientation_t {forward_forward, forward_reverse, reverse_forward};
enum class sam_results_t {undefined, mapped, unmapped};
enum class mapping_type_t {lev = 0, indel1 = 1, indel2 = 2, indel_clipping = 3, clipping_indel = 4, mismatches_clipping = 5, 
	clipping_mismatches = 6, clipping_clipping = 7, none = 8};

inline genome_t revert_direction(genome_t dir) { return dir == genome_t::direct ? genome_t::rev_comp : genome_t::direct; }

inline genome_t get_direction_using_orientation(genome_t dir, mapping_orientation_t mapping_orientation) {
	if (mapping_orientation == mapping_orientation_t::forward_forward) {
		return dir;
	}
	else if (mapping_orientation == mapping_orientation_t::forward_reverse) {
		return (dir == genome_t::direct) ? genome_t::rev_comp : genome_t::direct;
	}
	else if (mapping_orientation == mapping_orientation_t::reverse_forward) {
		// for now the same as co forward-reverse
		return (dir == genome_t::direct) ? genome_t::rev_comp : genome_t::direct;
	}
	else
	{
		// for now the same as co forward-reverse
		return (dir == genome_t::direct) ? genome_t::rev_comp : genome_t::direct;
	}
}

const uint32_t error_unknown = 255;

const uint32_t lev_alg_thr = 7;

const uint32_t shortest_min_read_len = 70;
const double largest_frac_errors = 0.05;
const uint32_t max_substr_len = 8;
const uint32_t sens_malicious_group_factor = 1;

typedef uint64_t read_id_t;
const read_id_t empty_read_id = ~((read_id_t) 0);

const uint32_t max_fastq_rec_length = 10240;

const string PRJ_CFG_FILE = "project.whisper";

// Index extenstions
const string EXT_DEF_IDX		 = ".whisper_idx";
const string EXT_REF_SEQ_DESC    = EXT_DEF_IDX + ".ref_seq_desc";
const string EXT_REF_SEQ_DIR     = EXT_DEF_IDX + ".ref_seq_dir";
const string EXT_REF_SEQ_RC      = EXT_DEF_IDX + ".ref_seq_rc";
const string EXT_REF_SEQ_DIR_PCK = EXT_DEF_IDX + ".ref_seq_dir_pck";
const string EXT_REF_SEQ_RC_PCK  = EXT_DEF_IDX + ".ref_seq_rc_pck";
const string EXT_SA_DIR          = EXT_DEF_IDX + ".sa_dir";
const string EXT_SA_RC           = EXT_DEF_IDX + ".sa_rc";
const string EXT_LUT_SHORT_DIR   = EXT_DEF_IDX + ".lut_short_dir";
const string EXT_LUT_SHORT_RC    = EXT_DEF_IDX + ".lut_short_rc";
const string EXT_LUT_LONG_DIR    = EXT_DEF_IDX + ".lut_long_dir";
const string EXT_LUT_LONG_RC     = EXT_DEF_IDX + ".lut_long_rc";
const string EXT_VCF			 = EXT_DEF_IDX + ".var";
const string EXT_REF_SNP		 = EXT_DEF_IDX + ".ref_snp";

// Project extensions
const string EXT_DEF_PRO		 = ".whisper_prj";
const string EXT_READS_BIN		 = EXT_DEF_PRO + ".reads_bin";
const string EXT_ID_STORE        = EXT_DEF_PRO + ".id_store";
const string EXT_BIN_STATS       = EXT_DEF_PRO + ".bin_stats";
const string EXT_MAPPING_GROUP   = EXT_DEF_PRO + ".mapping_group";


// Index markers
const string MARKER_DEF_IDX         = "MAP1";
const string MARKER_REF_SEQ_DESC    = MARKER_DEF_IDX + "RSDE"; 
const string MARKER_REF_SEQ_DIR     = MARKER_DEF_IDX + "RSDI"; 
const string MARKER_REF_SEQ_RC      = MARKER_DEF_IDX + "RSRC"; 
const string MARKER_REF_SEQ_DIR_PCK = MARKER_DEF_IDX + "RSDI_PCK"; 
const string MARKER_REF_SEQ_RC_PCK  = MARKER_DEF_IDX + "RSRC_PCK"; 
const string MARKER_SA_DIR          = MARKER_DEF_IDX + "SA_DIR"; 
const string MARKER_SA_RC           = MARKER_DEF_IDX + "SA_RC"; 
const string MARKER_LUT_SHORT_DIR   = MARKER_DEF_IDX + "LUT_SHORT_DIR"; 
const string MARKER_LUT_SHORT_RC    = MARKER_DEF_IDX + "LUT_SHORT_RC"; 
const string MARKER_LUT_LONG_DIR    = MARKER_DEF_IDX + "LUT_LONG_DIR"; 
const string MARKER_LUT_LONG_RC     = MARKER_DEF_IDX + "LUT_LONG_RC"; 
const string MARKER_VCF             = MARKER_DEF_IDX + "VAR";
const string MARKER_REF_SNP			= MARKER_DEF_IDX + "REF_SNP";

// Project markers
const string MARKER_DEF_PRO         = "MAP1_PRO";
const string MARKER_READS_BIN       = MARKER_DEF_PRO + "READS_BIN";
const string MARKER_ID_STORE        = MARKER_DEF_PRO + "ID_STORE";
const string MARKER_BIN_STATS       = MARKER_DEF_PRO + "BIN_STATS";
const string MARKER_MAPPING_GROUP   = MARKER_DEF_PRO + "MAPPING_GROUP";


// Running stats
const uint32_t STAT_READ_SPLITTING	  = 0x1000;

const uint32_t STAT_PRELOADING		  = 0x1800;

const uint32_t STAT_SORTING_BASE	  = 0x2000;
const uint32_t STAT_SORTING_TOTAL	  = STAT_SORTING_BASE + 0xfff;

const uint32_t STAT_BINS_WRITER_BASE  = 0x3000;
const uint32_t STAT_BINS_WRITER_TOTAL = STAT_BINS_WRITER_BASE + 0xfff;

const uint32_t STAT_MAPPING_GROUP_WRITER_BASE  = 0x4000;
const uint32_t STAT_MAPPING_GROUP_WRITER_TOTAL = STAT_MAPPING_GROUP_WRITER_BASE + 0xfff;

const uint32_t STAT_TIME_BASE	  = 0x5000;
const uint32_t STAT_TIME_TOTAL	  = STAT_TIME_BASE + 0xfff;
const uint32_t STAT_TIME_SPLIT	  = STAT_TIME_BASE + 0xffe;
const uint32_t STAT_TIME_MAIN	  = STAT_TIME_BASE + 0xffd;
const uint32_t STAT_TIME_POST     = STAT_TIME_BASE + 0xffc;

const uint32_t STAT_DIR_EXACT_MATCH	          = 0x6000;
const uint32_t STAT_RC_EXACT_MATCH			  = 0x6010;
const uint32_t STAT_DIR_OR_RC_EXACT_MATCH	  = 0x6020;
const uint32_t STAT_DIR_AND_RC_EXACT_MATCH	  = 0x6030;

const uint32_t STAT_MISMATCHES_TO_DIR_BASE	  = 0x7000;
const uint32_t STAT_MISMATCHES_TO_RC_BASE	  = 0x8000;

//const uint32_t STAT_DIR_ONE_MISMATCH		  = 0X9010;
//const uint32_t STAT_RC_ONE_MISMATCH			  = 0X9020;
const uint32_t STAT_SELECTED_READS_BASE		  = 0X9000;
const uint32_t STAT_CF_ACCEPTED				  = 0XA000;
const uint32_t STAT_CF_DISCARDED			  = 0XB000;
const uint32_t STAT_CF_ACCEPTED_IN_MALICIOUS_GROUPS	  = 0XC000;
const uint32_t STAT_CF_DISCARDED_IN_MALICIOUS_GROUPS  = 0XD000;
const uint32_t STAT_LEV_POSITIVE			   = 0XE000;
//const uint32_t STAT_LEV_NEGATIVE			   = 0XF000;
const uint32_t STAT_INDEL_MAPPING_FOUND	       = 0XF000;
const uint32_t STAT_INDEL_MAPPING_TESTS		   = 0X10000;

const uint32_t STAT_TIME_THR_RES_WRITER        = 0x20000;
const uint32_t STAT_TIME_THR_FASTQ_READER      = 0x20001;
const uint32_t STAT_TIME_THR_SPLITTER          = 0x20002;
const uint32_t STAT_TIME_THR_BIN_WRITER_BASE   = 0x21000;
const uint32_t STAT_TIME_THR_BIN_READER_BASE   = 0x22000;
const uint32_t STAT_TIME_THR_MAPPING_CORE_BASE = 0x23000;
const uint32_t STAT_TIME_THR_FASTQ_READER_PP   = 0x24000;
const uint32_t STAT_TIME_THR_RES_READER        = 0x24001;
const uint32_t STAT_TIME_THR_SAM_GENERATOR     = 0x24002;
const uint32_t STAT_TIME_THR_SAM_SORTING       = 0x24003;
const uint32_t STAT_TIME_THR_SAM_PROCESSING    = 0x24004;

const uint32_t STAT_PUSH_RESULTS               = 0x25000;
const uint32_t STAT_PUSH_RESULTS_UNQ           = 0x25001;
const uint32_t STAT_PUSH_RESULTS_SEND_BYTES    = 0x25002;
const uint32_t STAT_PUSH_RESULTS_ADDED_BYTES   = 0x25003;

const uint32_t STAT_MAPPED_SE_SINGLE		   = 0x26000;
const uint32_t STAT_MAPPED_SE_NONE			   = 0x26001;

const uint32_t STAT_MAPPED_PE_PAIR_UNQ		   = 0x26002;
const uint32_t STAT_MAPPED_PE_PAIR_MORE		   = 0x26003;

const uint32_t STAT_MAPPED_PE_ERRORS		   = 0x26006;

const uint32_t STAT_MAPPED_PE_INDEPENDENT_BOTH = 0x2600A;
const uint32_t STAT_MAPPED_PE_INDEPENDENT_SINGLE = 0x2600B;
const uint32_t STAT_MAPPED_PE_INDEPENDENT_NONE = 0x2600C;
const uint32_t STAT_MAPPED_PE_INDEPENDENT_ERRORS = 0x2600D;

const uint32_t STAT_MAPPED_PE_METHOD = 0x26010;

const uint32_t STAT_MAPPED_PE_HISTO = 0x27000;
const uint32_t STAT_MAPPED_PE_HISTO_CLIPPED = 0x28000;

const uint32_t STAT_SHORT_INDEL_REFINEMENTS_LEV = 0x28500;
const uint32_t STAT_SHORT_INDEL_REFINEMENTS_MAPPING = 0x28501;

const uint32_t STAT_MAPPED_PE_HISTO_VAR_INS_LONG = 0x29000;
const uint32_t STAT_MAPPED_PE_HISTO_VAR_DEL_LONG = 0x29500;
const uint32_t STAT_MAX_LEN_INDELS = 100;

const uint32_t STAT_BIN_TIMES = 0x30000;

const uint32_t STAT_READS_LEN = 0x31000;
const uint32_t STAT_READS_NS = 0x31300;
const uint32_t STAT_READS_LEN_WO_NS = 0x31600;

const uint32_t STAT_TIME_THR_PP_PARTS = 0x40000;

const uint32_t STAT_MAX_READ_LEN = 512;

#define MIN(_a, _b)			(((_a) < (_b)) ? (_a) : (_b))	
#define MAX(_a, _b)			(((_a) > (_b)) ? (_a) : (_b))	
#define MIN3(_a, _b, _c)	(MIN(_a, MIN(_b, _c)))	
#define MAX3(_a, _b, _c)	(MAX(_a, MAX(_b, _c)))	
#define DIST(_a, _b)		(((_a) > (_b)) ? ((_a) - (_b)) : ((_b) - (_a)))

#ifdef _DEBUG
#define A_memcpy	memcpy
#define A_memset	memset
#endif

#endif

// EOF
