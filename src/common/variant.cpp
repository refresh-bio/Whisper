// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 2.0
// Date    : 2018-10-18
// License : GNU GPL 3
// *******************************************************************************************

#include "../common/defs.h"
#include "../common/types.h"
#include "../common/variant.h"
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>

#ifdef ENABLE_VCF_VARIANTS

// *******************************************************************************************
CVariantDB::CVariantDB()
{
	max_lut_pos = 0;
	max_pos = 0;
}

// *******************************************************************************************
CVariantDB::~CVariantDB()
{

}

// *******************************************************************************************
void CVariantDB::Clear()
{

}

// *******************************************************************************************
bool CVariantDB::PrepareLUTs(uint32_t _max_pos)
{
	max_pos = _max_pos;
	max_lut_pos = (max_pos + 1) / lut_block_size;

	v_lut.clear();
	v_lut.resize(max_lut_pos+1ull, lut_type_t{ (uint32_t) v_snp.size()+1, (uint32_t)v_ins.size() + 1, (uint32_t)v_del.size() + 1, (uint32_t)v_sv.size() + 1});

	// Insert SNP into LUT
	for (uint32_t i = 0; i < v_snp.size(); ++i)
	{
		uint32_t cur_lut_pos = v_snp[i].pos / lut_block_size;
		v_lut[cur_lut_pos].snp = min(v_lut[cur_lut_pos].snp, cur_lut_pos);
	}

	for (uint32_t i = 0; i < v_ins.size(); ++i)
	{
		uint32_t cur_lut_pos = v_ins[i].pos / lut_block_size;
		v_lut[cur_lut_pos].ins = min(v_lut[cur_lut_pos].ins, cur_lut_pos);
	}

	for (uint32_t i = 0; i < v_del.size(); ++i)
	{
		uint32_t cur_lut_pos = v_del[i].pos / lut_block_size;
		v_lut[cur_lut_pos].del = min(v_lut[cur_lut_pos].del, cur_lut_pos);
	}

	for (uint32_t i = 0; i < v_sv.size(); ++i)
	{
		uint32_t cur_lut_pos = v_snp[i].pos / lut_block_size;
		v_lut[cur_lut_pos].sv = min(v_lut[cur_lut_pos].sv, cur_lut_pos);
	}

	// Fill empty spaces
	for (int i = (int) v_lut.size() - 2; i >= 0; --i)
	{
		v_lut[i].snp = min(v_lut[i].snp, v_lut[i + 1ull].snp);
		v_lut[i].ins = min(v_lut[i].ins, v_lut[i + 1ull].ins);
		v_lut[i].del = min(v_lut[i].del, v_lut[i + 1ull].del);
		v_lut[i].sv = min(v_lut[i].sv, v_lut[i + 1ull].sv);
	}

	return true;
}

// *******************************************************************************************
bool CVariantDB::InsertSNP(bool is_ref, uint8_t symbol, float freq, ref_pos_t pos)
{
	v_snp.push_back(variant_snp_t(
		(uint8_t)is_ref,
		DNA2Bin(symbol),
		(uint16_t) (freq * 65535),
		pos));

	return false;
}

// *******************************************************************************************
bool CVariantDB::InsertIns(const string &symbols, float freq, ref_pos_t pos)
{
	uint32_t data_pos = (uint32_t) v_ins_data.size();

	for (auto c : symbols)
		v_ins_data.push_back(DNA2Bin(c));

	v_ins.push_back(variant_ins_t(data_pos, (uint32_t) symbols.size(), (uint16_t)(freq * 65535), (uint32_t) pos));

	return true;
}

// *******************************************************************************************
bool CVariantDB::InsertDel(uint32_t len, float freq, ref_pos_t pos)
{
	v_del.push_back(variant_del_t(len, (uint16_t)(freq * 65535), pos));

	return true;
}

// *******************************************************************************************
bool CVariantDB::InsertSV(const string &symbols, uint32_t del_len, float freq, ref_pos_t pos)
{
	uint32_t data_pos = (uint32_t) v_sv_data.size();

	for (auto c : symbols)
		v_ins_data.push_back(DNA2Bin(c));

	v_sv.push_back(variant_sv_t(data_pos, (uint32_t) symbols.size(), del_len, (uint16_t)(freq * 65535), (uint32_t) pos));

	return true;
}

// *******************************************************************************************
uint8_t* CVariantDB::GetInsDesc(uint32_t ins_desc_pos)
{
	return v_ins_data.data() + ins_desc_pos;
}

// *******************************************************************************************
void CVariantDB::sort_db()
{
	sort(v_snp.begin(), v_snp.end(), [](variant_snp_t x, variant_snp_t y) {
		if (x.pos != y.pos)
			return x.pos < y.pos;
		if (x.is_ref != y.is_ref)
			return (bool)x.is_ref;
		return x.symbol < y.symbol;
	}
	);

	sort(v_ins.begin(), v_ins.end(), [](variant_ins_t x, variant_ins_t y) {
		return x.pos < y.pos;
	}
	);

	sort(v_del.begin(), v_del.end(), [](variant_del_t x, variant_del_t y) {
		return x.pos < y.pos;
	}
	);

	sort(v_sv.begin(), v_sv.end(), [](variant_sv_t x, variant_sv_t y) {
		return x.pos < y.pos;
	}
	);
}

// *******************************************************************************************
void CVariantDB::prepare_snp_ref(uchar_t *ref_seq_data, uint32_t ref_seq_size)
{
	v_snp_ref.resize(ref_seq_size);


	for (uint32_t i = 0; i < ref_seq_size * 2; ++i)
	{
		uchar_t ref_symbol = GET_SYMBOL(ref_seq_data, i);
		uint8_t code = 0;

		if (ref_symbol != sym_code_N_ref)
			code = (uint8_t)(1u << ref_symbol);

		v_snp_ref[i / 2] = code << (4 * ((i + 1) % 2));
	}

	for (auto x : v_snp)
	{
		uint8_t add_code = (uint8_t)(1u << x.symbol);
		v_snp_ref[x.pos / 2] |= add_code << (4 * ((x.pos + 1) % 2));
	}
}

// *******************************************************************************************
bool CVariantDB::Finalize(uchar_t *ref_seq_data, uint32_t ref_seq_size)
{
	sort_db();
	prepare_snp_ref(ref_seq_data, ref_seq_size);

	return true;
}

// *******************************************************************************************
bool CVariantDB::Serialize(shared_ptr<CMapperFile> f_vcf, shared_ptr<CMapperFile> f_snp_ref)
{
	// *** Variant databases 
	// Sizes
	f_vcf->WriteUint32((uint32_t)v_snp.size());
	f_vcf->WriteUint32((uint32_t)v_ins.size());
	f_vcf->WriteUint32((uint32_t)v_ins_data.size());
	f_vcf->WriteUint32((uint32_t)v_del.size());
	f_vcf->WriteUint32((uint32_t)v_sv.size());
	f_vcf->WriteUint32((uint32_t)v_sv_data.size());

	// Data
	for (auto x : v_snp)
		f_vcf->Write(&x, sizeof(x));

	for(auto x : v_ins)
		f_vcf->Write(&x, sizeof(x));
	f_vcf->Write(v_ins_data.data(), v_ins_data.size());

	for (auto x : v_del)
		f_vcf->Write(&x, sizeof(x));

	for (auto x : v_sv)
		f_vcf->Write(&x, sizeof(x));
	f_vcf->Write(v_sv_data.data(), v_sv_data.size());

	// *** Reference sequence with SNPs
	f_snp_ref->WriteUint32((uint32_t)v_snp_ref.size());
	f_snp_ref->Write(v_snp_ref.data(), v_snp_ref.size());

	return true;
}

// *******************************************************************************************
bool CVariantDB::Deserialize(shared_ptr<CMapperFile> f_vcf, shared_ptr<CMapperFile> f_snp_ref)
{
	// *** Variant databases
	v_snp.clear();
	v_ins.clear();
	v_ins_data.clear();
	v_del.clear();
	v_sv.clear();
	v_sv_data.clear();

	// Load sizes
	uint32_t tmp;
	if (f_vcf->ReadUint32(tmp) < 0)		return false;
	v_snp.resize(tmp);

	if (f_vcf->ReadUint32(tmp) < 0)		return false;
	v_ins.resize(tmp);

	if (f_vcf->ReadUint32(tmp) < 0)		return false;
	v_ins_data.resize(tmp);

	if (f_vcf->ReadUint32(tmp) < 0)		return false;
	v_del.resize(tmp);

	if (f_vcf->ReadUint32(tmp) < 0)		return false;
	v_sv.resize(tmp);

	if (f_vcf->ReadUint32(tmp) < 0)		return false;
	v_sv_data.resize(tmp);
	
	// Load data
	for (auto &p : v_snp)
		f_vcf->Read(&p, sizeof(p));

	for (auto &p : v_ins)
		f_vcf->Read(&p, sizeof(p));

	f_vcf->Read(v_ins_data.data(), v_ins_data.size());

	for (auto &p : v_del)
		f_vcf->Read(&p, sizeof(p));
	
	for (auto &p : v_sv)
		f_vcf->Read(&p, sizeof(p));

	f_vcf->Read(v_sv_data.data(), v_sv_data.size());

	// *** Reference sequence with SNPs
	if (f_snp_ref->ReadUint32(tmp) < 0)		return false;
	v_snp_ref.resize(tmp);
	
	f_snp_ref->Read(v_snp_ref.data(), tmp);

	PrepareLUTs(4000000000u);	// !!! Tymczasowo wart 4e9

	return true;
}

// *******************************************************************************************
// !!! Wstepna implementacja. Do sprawdzenia czy to dziala
pair<vector<variant_snp_t>::iterator, vector<variant_snp_t>::iterator> CVariantDB::FindRangeSNP(ref_pos_t first, ref_pos_t last)
{
/*	auto i_first = lower_bound(v_snp.begin() + v_lut[first / lut_block_size].snp, v_snp.begin() + v_lut[first / lut_block_size + 1ull].snp, first, 
		[](variant_snp_t x, ref_pos_t y) {
		return x.pos < y;
	});

	auto i_last = upper_bound(v_snp.begin() + v_lut[last / lut_block_size].snp, v_snp.begin() + v_lut[last / lut_block_size + 1ull].snp, last, 
		[](ref_pos_t x, variant_snp_t y) {
		return x < y.pos;
	});
*/
	auto i_first = lower_bound(v_snp.begin(), v_snp.end(), first, [](variant_snp_t x, ref_pos_t y) {
		return x.pos < y;
		});

	auto i_last = upper_bound(v_snp.begin(), v_snp.end(), last, [](ref_pos_t x, variant_snp_t y) {
		return x < y.pos;
		});

	return make_pair(i_first, i_last);
}

// *******************************************************************************************
// !!! Wstepna implementacja. Do sprawdzenia czy to dziala
pair<vector<variant_ins_t>::iterator, vector<variant_ins_t>::iterator> CVariantDB::FindRangeIns(ref_pos_t first, ref_pos_t last)
{
	auto i_first = lower_bound(v_ins.begin(), v_ins.end(), first, [](variant_ins_t x, ref_pos_t y) {
		return x.pos < y;
	});

	auto i_last = upper_bound(v_ins.begin(), v_ins.end(), last, [](ref_pos_t x, variant_ins_t y) {
		return x < y.pos;
	});

	return make_pair(i_first, i_last);
}

// *******************************************************************************************
// !!! Wstepna implementacja. Do sprawdzenia czy to dziala
pair<vector<variant_del_t>::iterator, vector<variant_del_t>::iterator> CVariantDB::FindRangeDel(ref_pos_t first, ref_pos_t last)
{
	auto i_first = lower_bound(v_del.begin(), v_del.end(), first, [](variant_del_t x, ref_pos_t y) {
		return x.pos < y;
	});

	auto i_last = upper_bound(v_del.begin(), v_del.end(), last, [](ref_pos_t x, variant_del_t y) {
		return x < y.pos;
	});

	return make_pair(i_first, i_last);
}

// *******************************************************************************************
// !!! Wstepna implementacja. Do sprawdzenia czy to dziala
pair<vector<variant_sv_t>::iterator, vector<variant_sv_t>::iterator> CVariantDB::FindRangeSV(ref_pos_t first, ref_pos_t last)
{
	auto i_first = lower_bound(v_sv.begin(), v_sv.end(), first, [](variant_sv_t x, ref_pos_t y) {
		return x.pos < y;
	});

	auto i_last = upper_bound(v_sv.begin(), v_sv.end(), last, [](ref_pos_t x, variant_sv_t y) {
		return x < y.pos;
	});

	return make_pair(i_first, i_last);
}

#endif

// EOF
