// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.0
// Date    : 2017-12-24
// License : GNU GPL 3
// *******************************************************************************************


#ifndef _SAM_H
#define _SAM_H

#include "../common/defs.h"
#include "../common/queue.h"
#include "../common/mmgr.h"
#include "../common/stats.h"
#include "../common/types.h"
#include "../common/params.h"
#include "reads.h"
#include "../common/joiner_mgr.h"
#include "soft_clipping.h"
#include "ref_desc.h"
#include "LevMyers.h"
#include "mapping_evaluation.h"
#include <random>

// ************************************************************************************
class CSamGenerator
{
	CParams *params;
	CReference *reference;

	CMemoryMonitor *mem_monitor;
	CRegisteringQueue<mapped_reads_t> *q_map_reads;
	CMemoryPool<uchar_t> *mp_fastq_blocks;
	CMemoryPool<uchar_t> *mp_sam_parts;
	CMemoryPool<uchar_t> *mp_ext_cigar;
	CMemoryPool<uchar_t> *mp_cigar;
	CMemoryPool<uchar_t> *mp_mdz;
	CRegisteringQueue<sam_block_t> *q_sam_blocks;

	CRunningStats *running_stats;

	CReadsReader *reads_reader;
	CReadsReader *reads_readers[2];
	CRefSeqDesc *ref_seq_desc;

	LevMyers* levMyers64;
	LevMyers* levMyers128;
	LevMyers* levMyers256;

	CSoftClipping* softClipping;

	InsertSizeModel* insertSizeModel;

	CProgress *progress;

	CGzipMember *gzip;
	
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

	uchar_t *data_srt;
	uchar_t *data_srt_end;
	uint64_t *ptrs;

	uint32_t max_mate_distance;
	double max_mate_edit_distance_frac;
	uchar_t *tmp_read_sequence;
	uchar_t *tmp_ref_sequence;

	uchar_t *mapped_part;
	uint64_t mapped_part_size;
	uint64_t mapped_part_pos;
	uint64_t mapped_part_reserve;
	
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

	std::vector<int64_t> stat_mapped_pe_histo;
	std::vector<int64_t> stat_mapped_pe_histo_clipped;
	std::vector<int64_t> stat_mapped_pe_method;

	vector<struct mapping_desc_t> mapping_desc[2];
	std::mt19937 engine;

	void process_group_se();
	void process_group_pe();

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

	void generate_unique_hits(std::vector<mapping_desc_t>& mapping_desc, std::vector<mapping_desc_t*>& hits);

	bool find_duplicate_pairs(const std::vector<mapping_pair_t>& mapping_pairs);

	template <MatchingMethod::Enumeration method>
	bool find_pairs(
		std::vector<mapping_desc_t> mapping_desc[], 
		std::vector<mapping_desc_t*> hits[],
		uchar_t *id[2], 
		uchar_t *sequence[2], 
		uint32_t sequence_len[2], 
		uchar_t *quality[2], 
		std::vector<mapping_pair_t>& mapping_pairs,
		std::pair<int, int>& editDistanceRange);

	bool find_single(
		std::vector<mapping_desc_t*> hits[], 
		uchar_t *id[2], 
		uchar_t *sequence[2], 
		uint32_t sequence_len[2], 
		uchar_t *quality[2]);

	void fill_cigar_with_lev(mapping_desc_t& mapping_desc, uchar_t *id, uchar_t *sequence, uint32_t sequence_len, uchar_t* quality, bool distant, bool set_affine_score);

	void extract_cigar_and_mdz(const uchar_t * ext_cigar, uchar_t *& cigar, uchar_t *& mdz);

	inline uint32_t trim_id(uchar_t *id, uint32_t id_len);
	inline void append_part(uchar_t *s, uint32_t s_len, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void append_part_rev(uchar_t *s, uint32_t s_len, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void append_part_rev_comp(uchar_t *s, uint32_t s_len, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void append_part(string s, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void append_part(uchar_t* s, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void append_part(int32_t x, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void append_part(double x, int prec, uchar_t *mapped_part, uint64_t &mapped_part_pos, bool add_tab);
	inline void close_line(uchar_t *mapped_part, uint64_t &mapped_part_pos);

	void push_mapped_part();
	
	inline void store_mapping_result(uchar_t *&buffer, uint64_t &pos,
		uchar_t *id, uint32_t id_len, uint32_t flag, string ref_seq_name, uint32_t ref_seq_pos, uint32_t mapping_quality, uchar_t* cigar,
		string mate_ref_seq_name, uint32_t mate_ref_seq_pos, int32_t template_len, genome_t dir, uchar_t *sequence, uint32_t sequence_len,
		uchar_t *qualities, uint32_t edit_distance, uchar_t* mdz, double score, MatchingMethod method);

	void convert_to_rev_comp(uchar_t *dest, uchar_t *src, uint32_t len);
	void copy_direct(uchar_t *dest, uchar_t *src, uint32_t len);
	void ref_copy_direct(uchar_t *dest, uchar_t *src, ref_pos_t pos, uint32_t len);

public:
	CSamGenerator(CParams *params, CObjects *objects, CRefSeqDesc *_ref_seq_desc);
	~CSamGenerator();

	void operator()();
};

#endif

// EOF
