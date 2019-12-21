#pragma once
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

#ifdef ENABLE_VCF_VARIANTS
#include "../common/defs.h"
#include "../common/types.h"
#include "../common/variant.h"
#include <string>
#include <vector>

using namespace std;

bool vcf_construct(const string index_name, const string dest_dir, const vector<string> &vcf_names);
bool load_ref_seq_desc(const string index_name, const string dest_dir, vector<struct seq_desc_t> &seq_desc);
	bool process_vcf_file(CVariantDB &variant_db, const vector<struct seq_desc_t> &seq_desc, const vector<string> &vcf_names,
		uchar_t *ref_seq_data, uint32_t ref_seq_size);

bool parse_vcf_record(FILE *in, string &chr, uint32_t &pos, string &id, string &ref, string &alt, string &qual, string &filter, string &info);
bool read_vcf_field(FILE *in, string &s, bool &was_eol);
bool move_to_eol(FILE *in);
void separate_string(string &s, vector<string> &v_tokens, char separator);
#endif

// EOF
