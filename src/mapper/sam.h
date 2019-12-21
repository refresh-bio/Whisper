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


#ifndef _SAM_H
#define _SAM_H

#include "../common/defs.h"
#include "queue.h"
#include "mmgr.h"
#include "stats.h"
#include "../common/types.h"
#include "params.h"
#include "reads.h"
#include "../common/variant.h"
#include "joiner_mgr.h"
#include "soft_clipping.h"
#include "indels.h"
#include "ref_desc.h"
#include "LevMyers.h"
#include "mapping_evaluation.h"
#include "sam_part.h"
#include "variant_match.h"

#include <random>
#include <functional>

//#define VALIDATE_CIGAR
//#define REFINE_SCORE_SNP

#ifdef COLLECT_STATS
#define MEASURE_PP_TIMES
#endif

// ************************************************************************************
class CSamGenerator
{
	const array<uchar_t, 2> a_AS{ 'A', 'S' };
	const array<uchar_t, 2> a_MD{ 'M', 'D' };
	const array<uchar_t, 2> a_NM{ 'N', 'M' };
	const array<uchar_t, 2> a_RG{ 'R', 'G' };
	const array<uchar_t, 2> a_XM{ 'X', 'M' };
	const uchar_t uc_star[2] = { '*', 0 };
	const int mapping_heap_size = 128;

	CParams *params;
	CReference *reference;
#ifdef ENABLE_VCF_VARIANTS
	CVariantDB *variant_db;
#endif

	FILE *f_out[2];

	CMemoryMonitor *mem_monitor;
	CRegisteringQueue<mapped_reads_t> *q_map_reads;
	CMemoryPool<uchar_t> *mp_fastq_blocks;
	CMemoryPool<uchar_t> *mp_ext_cigar;
	CMemoryPool<uchar_t> *mp_cigar;
	CMemoryPool<uint32_t> *mp_cigar_bin;
	CMemoryPool<uchar_t> *mp_mdz;
	CRegisteringQueue<sam_block_t> *q_sam_blocks;
	CPtrPool *ptr_pool;

	array<int64_t, STAT_MAX_READ_LEN + 1> stat_read_len;
	array<int64_t, STAT_MAX_READ_LEN + 1> stat_read_len_wo_Ns;
	array<int64_t, STAT_MAX_READ_LEN + 1> stat_read_Ns;
	CRunningStats *running_stats;

#ifdef MEASURE_PP_TIMES
	array<double, 16> part_times{ 0 };
#endif

	CReadsReader *reads_reader;
	CReadsReader *reads_readers[2];
	CRefSeqDesc *ref_seq_desc;

	LevMyers* levMyers64;
	LevMyers* levMyers128;
	LevMyers* levMyers256;

	CSoftClipping* softClipping;
	CIndelMatching* indelMatching;
#ifdef ENABLE_VCF_VARIANTS
	CVariantMatch* variantMatch;
#endif
	InsertSizeModel* insertSizeModel;

	CProgress *progress;

	uint32_t verbosity_level;
	uint32_t id_bits_subgroup;
	uint32_t id_bits_local;

	mapping_mode_t mapping_mode;
	mapping_orientation_t mapping_orientation;
	mapped_reads_t mapped_reads;
	uint32_t max_no_mappings;
	uint32_t mapping_counter_size;
	uint32_t id_bytes;

	uint32_t gzipped_SAM_level;
	bool store_BAM;

	uchar_t *data_srt;
	uchar_t *data_srt_end;
	vector<pair<uint64_t, uint64_t>> ptrs;

	uint32_t max_mate_distance;
	double max_mate_edit_distance_frac;
	uchar_t *tmp_read_sequence;
	uchar_t *tmp_ref_sequence;

	CSamPart *sam_part;
	CBamPart *bam_part;
	
	uint64_t stored_mapped;
	uint64_t stored_unmapped;
	int64_t read_bytes;

	int64_t stat_mapped_se_single;
	int64_t stat_mapped_se_none;
	int64_t stat_mapped_pe_pair_unq;
	int64_t stat_mapped_pe_pair_more;

	int64_t stat_mapped_pe_independent_both;
	int64_t stat_mapped_pe_independent_single;
	int64_t stat_mapped_pe_independent_none;
	int64_t stat_mapped_pe_independent_errors;
	int64_t stat_mapped_pe_errors;
	int64_t stat_refined_short_indels_lev;
	int64_t stat_refined_short_indels_mapping;

	std::vector<int64_t> stat_mapped_pe_histo;
	std::vector<int64_t> stat_mapped_pe_histo_clipped;
	std::vector<int64_t> stat_mapped_pe_method;

	vector<struct mapping_desc_t> mapping_desc[2];
	std::mt19937 engine;

	const function<bool(mapping_desc_t &, mapping_desc_t &)> mapping_desc_comparator = [](mapping_desc_t &x, mapping_desc_t &y) {
		if (x.err_edit_distance != y.err_edit_distance)
			return x.err_edit_distance < y.err_edit_distance;

		// Use hash instead of plain raw_pos to pick mappings from "random" positions instead of prefering first chromosomes
		return (x.mapping.pos * 0x678dde6fu) < (y.mapping.pos * 0x678dde6fu);
	};
	
	void process_group_se();
	void process_group_pe(uint64_t group_id);

	void sort_mapping_results();
	void free_reads();

	void check_buffers_occupation();

	bool find_id(read_id_t read_id, uchar_t *data_ptr, uchar_t *data_srt_end);
	bool find_id_pair(read_id_t read_id, uchar_t *data_ptr, uchar_t *data_srt_end);
	bool adjust_pos_in_results(read_id_t read_id, uchar_t *&data_ptr);

	void store_mapped_read(uchar_t *id, uchar_t *sequence, uchar_t *plus, uchar_t *quality, 
		uint32_t id_len, uint32_t sequence_len, uint32_t plus_len, uint32_t quality_len,
		uchar_t *&data_ptr);
	void store_mapped_pair_reads(read_id_t read_id, uchar_t *id[2], uchar_t *sequence[2], uchar_t *plus[2], uchar_t *quality[2], 
		uint32_t id_len[2], uint32_t sequence_len[2], uint32_t plus_len[2], uint32_t quality_len[2],
		uchar_t *&data_ptr);

	void store_unmapped_read(uchar_t *id, uchar_t *sequence, uchar_t *plus, uchar_t *quality, 
		uint32_t id_len, uint32_t sequence_len, uint32_t plus_len, uint32_t quality_len,
		uchar_t *&data_ptr);


	bool update_model(read_id_t read_id, uchar_t *id[2], uchar_t *sequence[2],
		uint32_t id_len[2], uint32_t sequence_len[2], uchar_t *&data_ptr);

	void generate_unique_hits(std::vector<mapping_desc_t>& mapping_desc, std::vector<mapping_desc_t*>& hits);

	bool find_duplicate_pairs(const std::vector<mapping_pair_t>& mapping_pairs);

	bool is_clipped(mapping_type_t t)
	{
		return
			t == mapping_type_t::clipping_clipping ||
			t == mapping_type_t::clipping_indel ||
			t == mapping_type_t::clipping_mismatches ||
			t == mapping_type_t::indel_clipping ||
			t == mapping_type_t::mismatches_clipping;
	}

	template <MatchingMethod::Enumeration method>
	bool find_pairs(
		std::vector<mapping_desc_t> mapping_desc[2], 
		std::vector<mapping_desc_t*> hits[2],
		uchar_t *id[2], 
		uchar_t *sequence[2], 
		uint32_t sequence_len[2], 
		uchar_t *quality[2], 
		std::vector<mapping_pair_t>& mapping_pairs,
		std::pair<int, int>& editDistanceRange);

	bool find_single(
		std::vector<mapping_desc_t*> hits[2], 
		uchar_t *id[2], 
		uchar_t *sequence[2], 
		uint32_t sequence_len[2], 
		uchar_t *quality[2]);

	void fill_cigar(mapping_desc_t& mapping_desc, uchar_t* id, uchar_t* sequence, uint32_t sequence_len, uchar_t* quality, bool distant, bool set_affine_score);
	void fill_cigar_with_lev(mapping_desc_t& mapping_desc, uchar_t *id, uchar_t *sequence, uint32_t sequence_len, uchar_t* quality, bool distant, bool set_affine_score);
	void fill_cigar_with_clipping(mapping_desc_t& mapping_desc, uchar_t* id, uchar_t* sequence, uint32_t sequence_len, uchar_t* quality, bool set_affine_score);
	void fill_cigar_with_indel_matching(mapping_desc_t& mapping_desc, uchar_t* id, uchar_t* sequence, uint32_t sequence_len, uchar_t* quality, bool set_affine_score);

	int extract_cigar_and_mdz(const uchar_t * ext_cigar, uchar_t * cigar, uchar_t * mdz);
	int extract_cigarbin_and_mdz(const uchar_t * ext_cigar, uint32_t * cigar_bin, uchar_t * mdz);

#ifdef VALIDATE_CIGAR
	bool validate_cigar(mapping_desc_t& md, uchar_t *id, uint32_t id_len, uchar_t* sequence, uint32_t sequence_len, genome_t direction);
#endif
#ifdef REFINE_SCORE_SNP
	void refine_score_snp(mapping_desc_t& md, uchar_t* id, uint32_t id_len, uchar_t* sequence, uint32_t sequence_len, genome_t direction);
#endif

	inline uint32_t trim_id(uchar_t *id, uint32_t id_len);

	void push_mapped_part();
	
	inline void store_mapping_result(
		uchar_t *id, uint32_t id_len, uint32_t flag, int32_t ref_seq_id, uint32_t ref_seq_pos, uint32_t mapping_quality, uchar_t* cigar, uint32_t *cigar_bin, uint32_t cigar_bin_len,
		int32_t mate_ref_seq_id, uint32_t mate_ref_seq_pos, int32_t template_len, genome_t dir, int32_t ref_length, uchar_t *sequence, uint32_t sequence_len,
		uchar_t *qualities, uint32_t edit_distance, uchar_t* mdz, double score, MatchingMethod method);

	void convert_to_rev_comp(uchar_t *dest, uchar_t *src, uint32_t len);
	void copy_direct(uchar_t *dest, uchar_t *src, uint32_t len);
	void ref_copy_direct(uchar_t *dest, uchar_t *src, ref_pos_t pos, uint32_t len);

	bool refine_ext_cigar(uchar_t* ext_cigar, uchar_t* ref_seq, uint32_t& edit_distance, uint32_t& num_events, double& score);

	constexpr int decode_indel(uchar_t x)
	{
		if (x > 128)
			return -(int) (x - 128);
		else
			return (int)x;
	}

	bool pushMappingHeap(std::vector<mapping_pair_t>& mapping_pairs, const mapping_pair_t& candidate, int heapSize);

public:
	CSamGenerator(CParams *params, CObjects *objects, CRefSeqDesc *_ref_seq_desc);
	~CSamGenerator();

	void operator()();
};

#endif

// EOF
