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
#include "vcf.h"
#include "../common/defs.h"
#include "../common/utils.h"
#include <algorithm>

// *******************************************************************************************
bool load_ref_seq_desc(const string index_name, const string dest_dir, vector<struct seq_desc_t> &seq_desc)
{
	seq_desc.clear();

	shared_ptr<CMapperFile> in_desc(new CMapperFile(EXT_REF_SEQ_DESC, MARKER_REF_SEQ_DESC, 1));
	if (!in_desc->OpenRead(dest_dir + index_name))
		return false;

	for (int32_t id = 0; !in_desc->Eof(); ++id)
	{
		struct seq_desc_t p;
		uint32_t seq_name_len;

		p.id = id;
		in_desc->Read(&seq_name_len, sizeof(uint32_t));
		char *s = new char[seq_name_len + 1ull];
		s[seq_name_len] = 0;
		in_desc->Read(s, seq_name_len);
		p.name = string(s);
		delete[] s;
		in_desc->Read(&(p.size), sizeof(uint64_t));
		in_desc->Read(&(p.pos_in_ref_seq), sizeof(uint64_t));
		in_desc->Read(&(p.no_initial_Ns), sizeof(uint64_t));
		in_desc->Read(&(p.no_final_Ns), sizeof(uint64_t));
		in_desc->Read(&(p.no_starting_Ns), sizeof(uint64_t));

		// Trim sequence name
		p.name.resize(find(p.name.begin(), p.name.end(), ' ') - p.name.begin());

		seq_desc.push_back(p);
	}

	return true;
}

// *******************************************************************************************
bool read_vcf_field(FILE *in, string &s, bool &was_eol)
{
	s.clear();
	was_eol = false;

	while (true)
	{
		int c = getc(in);
		if (c == EOF)
			return false;
		if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
		{
			if (s.empty())
				continue;
			else
			{
				was_eol = (c == '\r' || c == '\n');

				return true;
			}
		}

		s.push_back(c);
	}

	return true;
}

// *******************************************************************************************
bool move_to_eol(FILE *in)
{
	while (true)
	{
		int c = getc(in);

		if (c == EOF || c == '\n' || c == '\r')
			return true;
	}

	return false;
}

// *******************************************************************************************
bool parse_vcf_record(FILE *in, string &chr, uint32_t &pos, string &id, string &ref, string &alt, string &qual, string &filter, string &info)
{
	if (feof(in))
		return false;

	string s_pos;
	bool was_eol;

	if (!read_vcf_field(in, chr, was_eol))		return false;

	if (chr[0] == '#')		// Comment line
	{
		if (!was_eol)
			move_to_eol(in);
		chr.clear();

		return true;
	}

	if (!read_vcf_field(in, s_pos, was_eol))	return false;
	if (!read_vcf_field(in, id, was_eol))		return false;
	if (!read_vcf_field(in, ref, was_eol))		return false;
	if (!read_vcf_field(in, alt, was_eol))		return false;
	if (!read_vcf_field(in, qual, was_eol))		return false;
	if (!read_vcf_field(in, filter, was_eol))	return false;
	if (!read_vcf_field(in, info, was_eol))		return false;

	if (!was_eol)	
		move_to_eol(in);

	pos = stoi(s_pos);
	
	return true;
}

// *******************************************************************************************
bool process_vcf_file(CVariantDB &variant_db, const vector<struct seq_desc_t> &seq_desc, const vector<string> &vcf_names,
	uchar_t *ref_seq_data, uint32_t ref_seq_size)
{
	uint32_t n_variant = 0;

	for (auto fn : vcf_names)
	{
		FILE *in = fopen(fn.c_str(), "rb");
		if (!in)
		{
			std::cerr << "Cannot open " << fn << endl;
			return false;
		}

		std::cerr << "\nProcessing: " << fn << endl;

		string chr;
		uint32_t pos;
		string id;
		string ref;
		string alt;
		string qual;
		string filter;
		string info;

		vector<string> v_alt, v_info, v_caf;
		vector<float> v_freq;

		while (parse_vcf_record(in, chr, pos, id, ref, alt, qual, filter, info))
		{
			if (chr.empty())	// Comment line
				continue;

			// !!! To dziala tylko na naszych danych. Koniecznie trzeba to jakos uogolnic, np. podawac ten doklejany prefiks
//			if (chr[0] != 'c')
//				chr = "chr" + chr;

			auto p_chr = find_if(seq_desc.begin(), seq_desc.end(), [&chr](const seq_desc_t x) {
				return x.name == chr;
				}
			);

			if (p_chr != seq_desc.end())
			{
				uint32_t global_pos = pos - 1;

				global_pos += (uint32_t) p_chr->pos_in_ref_seq;
				global_pos += (uint32_t) p_chr->no_starting_Ns;
				global_pos -= (uint32_t) p_chr->no_initial_Ns;

				if (ref.size() == 1)
				{
					// Check reference symbol in ref. seq.
					uchar_t ref_symbol = GET_SYMBOL(ref_seq_data, global_pos);

					if (ref_symbol != DNA2Bin(ref[0]))
						cerr << "Reference symbol from VCF file does not match symbol from reference: " << chr << " " << pos << endl;

					separate_string(alt, v_alt, ',');
					separate_string(info, v_info, ';');
					v_caf.clear();

					bool is_ref_freq = false;

					for (auto x : v_info)
						if (x.substr(0, 4) == "CAF=" || x.substr(0, 3) == "AF=")
						{
							is_ref_freq = x[0] == 'C';
							string tmp = x.substr(3ll + (int) is_ref_freq, x.size() - (3ll + is_ref_freq));
							separate_string(tmp, v_caf, ',');
							break;
						}

					for (auto &x : v_caf)
						if (x == ".")
							x = "0.0";

					if (v_caf.size() == v_alt.size())
					{
						// Add frequency of the reference value
						double ref_freq = 1.0;
						for (auto &x : v_caf)
							ref_freq -= stof(x);
						variant_db.InsertSNP(true, ref[0], (float) ref_freq, global_pos);
					}
					else
						variant_db.InsertSNP(true, ref[0], !v_caf.empty() ? stof(v_caf[0]) : 0, global_pos);

					for (int i = 0; i < v_alt.size(); ++i)
					{
						string s_caf = (i + (int)is_ref_freq < (int)v_caf.size()) ? v_caf[i + (int64_t)is_ref_freq] : "0";

						if (v_alt[i].size() == 1)		// SNP
							variant_db.InsertSNP(false, v_alt[i][0], stof(s_caf), global_pos);
						else // insertion
						{
							uint32_t len = (uint32_t) v_alt[i].size();
							variant_db.InsertIns(v_alt[i], stof(s_caf), global_pos);
						}
					}
				}
				else
				{
					separate_string(alt, v_alt, ',');
					separate_string(info, v_info, ';');

					bool is_ref_freq = false;

					v_caf.clear();

					for (auto x : v_info)
						if (x.substr(0, 4) == "CAF=" || x.substr(0, 3) == "AF=")
						{
							is_ref_freq = x[0] == 'C';
							string tmp = x.substr(3ll + (int) is_ref_freq, x.size() - (3ll + (int) is_ref_freq));
							separate_string(tmp, v_caf, ',');
							break;
						}

					for (auto &x : v_caf)
						if (x == ".")
							x = "0.0";

					for (int i = 0; i < v_alt.size(); ++i)
					{
						string s_caf = (i + (int)is_ref_freq < (int)v_caf.size()) ? v_caf[i + (int64_t)is_ref_freq] : "0";

						if (v_alt[i].size() == 1)
							variant_db.InsertDel((uint32_t) ref.size(), stof(s_caf), global_pos);
						else
						{
							// Try simplify the variant, e.g.
							// TAAA -> AAAA can be stored as SNP
							// TAAA -> TA can be stored as AAA -> A at next position
							uint32_t tmp_pos = global_pos;
							auto tmp_ref = ref;
							auto tmp_alt = v_alt[i];

							// Try remove trailing symbols
							while (tmp_ref.size() > 1 && tmp_alt.size() > 1)
							{
								if (tmp_ref.back() == tmp_alt.back())
								{
									tmp_ref.pop_back();
									tmp_alt.pop_back();
								}
								else
									break;
							}

							//Try remove front symbols
							while (tmp_ref.size() > 1 && tmp_alt.size() > 1)
							{
								if (tmp_ref.front() == tmp_alt.front())
								{
									tmp_ref.erase(tmp_ref.begin());
									tmp_alt.erase(tmp_alt.begin());
									tmp_pos++;
								}
								else 
									break;
							}

							if (tmp_ref.size() == 1 && tmp_alt.size() == 1)
								variant_db.InsertSNP(false, tmp_alt[0], stof(s_caf), tmp_pos);
							else if(tmp_ref.size() == 1 && tmp_alt.size() > 1)
								variant_db.InsertIns(tmp_alt, stof(s_caf), tmp_pos);
							else if(tmp_ref.size() > 1 && tmp_alt.size() == 1)
								variant_db.InsertDel((uint32_t) ref.size(), stof(s_caf), tmp_pos);
							else
								variant_db.InsertSV(tmp_alt, (uint32_t) tmp_ref.size(), stof(s_caf), tmp_pos);
						}
					}
				}
			}

			if (++n_variant % 10000 == 0)
//			if (++n_variant % 10000 == 0 || n_variant > 0)
				std::cerr << n_variant << " -  "
				"  SNP: " << variant_db.v_snp.size() <<
				"  INS: " << variant_db.v_ins.size() <<
				"  INS data: " << variant_db.v_ins_data.size() <<
				"  DEL: " << variant_db.v_del.size() <<
				"  SV: " << variant_db.v_sv.size() <<
				"  SV data: " << variant_db.v_sv_data.size() <<
				"  chr: " << chr << 
				"  pos: " << pos <<
				"\r";
		}

		fclose(in);

		std::cerr << n_variant << " -  " <<
			"  SNP: " << variant_db.v_snp.size() <<
			"  INS: " << variant_db.v_ins.size() <<
			"  INS data: " << variant_db.v_ins_data.size() <<
			"  DEL: " << variant_db.v_del.size() <<
			"  SV: " << variant_db.v_sv.size() <<
			"  SV data: " << variant_db.v_sv_data.size() <<
			"\n";
	}

	variant_db.Finalize(ref_seq_data, ref_seq_size);

	std::cerr << "End of processing\n";

	return true;
}

// *******************************************************************************************
bool vcf_construct(const string index_name, const string dest_dir, const vector<string> &vcf_names)
{
	std::cerr << "\n***** Stage 9. SNP and short indel construction\n";

	// Load ref. seq. desc
	vector<struct seq_desc_t> seq_desc;

	if (!load_ref_seq_desc(index_name, dest_dir, seq_desc))
		return false;

	shared_ptr<CMapperFile> ref_seq(new CMapperFile(EXT_REF_SEQ_DIR_PCK, MARKER_REF_SEQ_DIR_PCK, 1));
	if (!ref_seq->OpenRead(dest_dir + index_name))
		return false;

	uint32_t ref_seq_size = (uint32_t) ref_seq->GetSize();
	uchar_t *ref_seq_data = new uchar_t[ref_seq_size];

	ref_seq->Read(ref_seq_data, ref_seq_size);

	CVariantDB variant_db;

	if (!process_vcf_file(variant_db, seq_desc, vcf_names, ref_seq_data, ref_seq_size))
		return false;

	shared_ptr<CMapperFile> out_desc(new CMapperFile(EXT_VCF, MARKER_VCF, 1));
	if (!out_desc->OpenWrite(dest_dir + index_name))
		return false;

	shared_ptr<CMapperFile> out_ref_snp(new CMapperFile(EXT_REF_SNP, MARKER_REF_SNP, 1));
	if (!out_ref_snp->OpenWrite(dest_dir + index_name))
		return false;

	variant_db.Serialize(out_desc, out_ref_snp);

	delete[] ref_seq_data;

	return true;
}

// *******************************************************************************************
void separate_string(string &s, vector<string> &v_tokens, char separator)
{
	v_tokens.clear();

	auto p = s.begin();

	while(true)
	{
		auto q = find(p + 1, s.end(), separator);
		v_tokens.push_back(string(p, q));
		
		if(q == s.end())
			break;
		else
			p = q + 1;
	}
}
#endif

// EOF
