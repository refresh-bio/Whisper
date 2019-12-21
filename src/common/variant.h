#pragma once
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

#include "../common/defs.h"
#include "../common/types.h"
#include <string>
#include <vector>

#ifdef ENABLE_VCF_VARIANTS
class CVariantDB
{
public:		// tymczasowo - do debugowania
	const uint32_t lut_block_size = 1024;
	uint32_t max_pos;
	uint32_t max_lut_pos;

	vector<variant_snp_t> v_snp;
	vector<variant_ins_t> v_ins;
	vector<variant_del_t> v_del;
	vector<variant_sv_t> v_sv;

	typedef struct {
		uint32_t snp;
		uint32_t ins;
		uint32_t del;
		uint32_t sv;
	} lut_type_t;

	vector<lut_type_t> v_lut;

	vector<uint8_t> v_ins_data;
	vector<uint8_t> v_sv_data;

	vector<uint8_t> v_snp_ref;

	void sort_db();
	void prepare_snp_ref(uchar_t *ref_seq_data, uint32_t ref_seq_size);

public: 
	CVariantDB();
	~CVariantDB();

	void Clear();
	bool InsertSNP(bool is_ref, uint8_t symbol, float freq, ref_pos_t pos);
	bool InsertIns(const string &symbols, float freq, ref_pos_t pos);
	bool InsertDel(uint32_t len, float freq, ref_pos_t pos);
	bool InsertSV(const string &symbols, uint32_t del_len, float freq, ref_pos_t pos);
	bool Finalize(uchar_t *ref_seq_data, uint32_t ref_seq_size);
	bool Serialize(shared_ptr<CMapperFile> f_vcf, shared_ptr<CMapperFile> f_snp_ref);
	bool Deserialize(shared_ptr<CMapperFile> f_vcf, shared_ptr<CMapperFile> f_snp_ref);

	bool PrepareLUTs(uint32_t _max_pos);

	pair<vector<variant_snp_t>::iterator, vector<variant_snp_t>::iterator> FindRangeSNP(ref_pos_t first, ref_pos_t last);
	pair<vector<variant_ins_t>::iterator, vector<variant_ins_t>::iterator> FindRangeIns(ref_pos_t first, ref_pos_t last);
	pair<vector<variant_del_t>::iterator, vector<variant_del_t>::iterator> FindRangeDel(ref_pos_t first, ref_pos_t last);
	pair<vector<variant_sv_t>::iterator, vector<variant_sv_t>::iterator> FindRangeSV(ref_pos_t first, ref_pos_t last);

	uint8_t* GetInsDesc(uint32_t ins_desc_pos);
};

#endif

// EOF
