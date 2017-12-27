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


#include "ref_seq.h"
#include "../common/defs.h"
#include "../common/utils.h"

#include <iostream>
#include <memory>

char sym_codes[256];

// ************************************************************************************
// Append a single file to the reference
bool ref_seq_append(vector<struct seq_desc_t> &seq_desc, char *seq, uint64_t &pos, const string file_name)
{
	FILE *in = fopen(file_name.c_str(), "rb");

	if(!in)
	{
		cout << "Cannot open file: " << file_name << "\n";
		return false;
	}

	cout << "Processing: " << file_name << " ";
	setvbuf(in, nullptr, _IOFBF, 1 << 24);

	while(!feof(in))
	{
		string name;
		uint64_t size = 0;
		uint64_t pos_in_ref_seq = pos;
		uint64_t no_initial_Ns = 0;
		uint64_t no_final_Ns = 0;
		uint64_t no_starting_Ns = SEPARATE_N_LEN;

		int c = getc(in);
		if(c != '>')
		{
			fclose(in);
			return false;
		}

		// Read sequence name
		while(!feof(in))
		{
			c = getc(in);
			if(c == '\n' || c == '\r')
				break;
			name.push_back(c);
		}

		// Read sequence data
		bool init_Ns_state = true;
		while(true)
		{
			c = getc(in);
			if(c == EOF)
				break;
			if(c == '>')
			{
				ungetc(c, in);
				break;
			}
			if(c == '\n' || c == '\r')
				continue;

			if((c == 'N' || c == 'n') && init_Ns_state)
			{
				++no_initial_Ns;
				continue;
			}

			if(init_Ns_state)
			{
				// Insert MAX_ERROR Ns
				for(int i = 0; i < SEPARATE_N_LEN; ++i)
					seq[pos++] = sym_code_N_ref;
				init_Ns_state = false;
			}

			seq[pos+size++] = sym_codes[c];
			if((size & ((1 << 24) - 1)) == 0)
				cout << ".";
		}

		while(size > 0)
		{
			if(seq[pos+size-1] == sym_code_N_ref)
			{
				--size;
				++no_final_Ns;
			}
			else
				break;
		}
		cout << "\n";

		if(size)
		{
			seq_desc.push_back(seq_desc_t(name, -1, size, pos_in_ref_seq, no_initial_Ns, no_final_Ns, no_starting_Ns));
			pos += size;
		}
	}

	fclose(in);

	return true;
}

// ************************************************************************************
// Construct reference sequence from the input files
bool ref_seq_dir_construct(const string index_name, const vector<string> &ref_seq_names, const string dest_dir, const string temp_dir)
{
	vector<struct seq_desc_t> seq_desc;

	cout << "\n***** Stage 1. Reference sequence construction (direct)\n";
	cout << "Processing " << ref_seq_names.size() << " file(s)\n";

	uint64_t estimated_size = 0;

	for (auto &fn : ref_seq_names)
	{
		FILE *f = fopen(fn.c_str(), "rb");
		if (f)
		{
			my_fseek(f, 0, SEEK_END);
			estimated_size += my_ftell(f);
			fclose(f);
		}
	}
	cout << "Total length of the input files: " << estimated_size << " bytes\n";

	estimated_size += 2 * BOUNDARY_N_LEN + (ref_seq_names.size() - 1) * SEPARATE_N_LEN;

	shared_ptr<char> ref_seq(new char[estimated_size]);
	uint64_t pos = 0;

	for(int i = 0; i < 256; ++i)
		sym_codes[i] = 15;
	sym_codes['A'] = sym_codes['a'] = sym_code_A;
	sym_codes['C'] = sym_codes['c'] = sym_code_C;
	sym_codes['G'] = sym_codes['g'] = sym_code_G;
	sym_codes['T'] = sym_codes['t'] = sym_code_T;
	sym_codes['N'] = sym_codes['n'] = sym_code_N_ref;

	// ***** Construct direct version of ref. seq.
	// Insert initial Ns
	for(int i = 0; i < BOUNDARY_N_LEN - SEPARATE_N_LEN; ++i)
		ref_seq.get()[pos++] = sym_code_N_ref;

	// Append sequences
	for(auto &fn : ref_seq_names)
		if(!ref_seq_append(seq_desc, ref_seq.get(), pos, fn))
			return false;

	// Insert final Ns
	for(int i = 0; i < BOUNDARY_N_LEN; ++i)
		ref_seq.get()[pos++] = sym_code_N_ref;

	// ***** Save files
	shared_ptr<CMapperFile> out_ref(new CMapperFile(EXT_REF_SEQ_DIR, MARKER_REF_SEQ_DIR, 1));
	if(!out_ref->OpenWrite(temp_dir + index_name))
		return false;

	shared_ptr<CMapperFile> out_desc (new CMapperFile(EXT_REF_SEQ_DESC, MARKER_REF_SEQ_DESC, 1));
	if(!out_desc->OpenWrite(dest_dir + index_name))
		return false;

	out_ref->Write(ref_seq.get(), pos);

	for(auto &p : seq_desc)
	{
		uint32_t seq_name_len = (uint32_t) p.name.size();
		out_desc->Write(&seq_name_len, sizeof(uint32_t));
		out_desc->Write(p.name.c_str(), p.name.size());
		out_desc->Write(&(p.size), sizeof(uint64_t));
		out_desc->Write(&(p.pos_in_ref_seq), sizeof(uint64_t));
		out_desc->Write(&(p.no_initial_Ns), sizeof(uint64_t));
		out_desc->Write(&(p.no_final_Ns), sizeof(uint64_t));
		out_desc->Write(&(p.no_starting_Ns), sizeof(uint64_t));
	}

	return true;
}

// ************************************************************************************
// Construct rev. comp. reference sequence relating on direct ref. seq.
bool ref_seq_rc_construct(const string index_name, const string temp_dir)
{
	cout << "\n***** Stage 2. Reference sequence construction (rev. comp.)\n";

	shared_ptr<CMapperFile> in(new CMapperFile(EXT_REF_SEQ_DIR, MARKER_REF_SEQ_DIR));

	if(!in->OpenRead(temp_dir + index_name))
	{
		cout << "Cannot open reference sequence (direct)\n";
		return false;
	}
	int64_t size = in->GetSize();

	shared_ptr<CMapperFile> out(new CMapperFile(EXT_REF_SEQ_RC, MARKER_REF_SEQ_RC));
	if(!out->OpenWrite(temp_dir + index_name))
	{
		cout << "Cannot create reference sequence (rev. comp.)\n";
		return false;
	}

	cout << "Reading...\n";
	shared_ptr<char> ref_seq(new char[size]);
	in->Read(ref_seq.get(), size);

	cout << "Transforming...\n";
	char *ptr = ref_seq.get();
	for(int64_t i = 0; i < (size+1)/2; ++i)
	{
		char b1 = ptr[i];
		char b2 = ptr[size-i-1];
		
		// complement (assumming encoding: A-0, C-1, G-2, T-3)
		if(b2 < sym_code_N_ref)
			ptr[i]        = 3 - b2;			
		else
			ptr[i]        = sym_code_N_ref;
		if(b1 < sym_code_N_ref)
			ptr[size-i-1] = 3 - b1;
		else
			ptr[size-i-1] = sym_code_N_ref;
	}

	cout << "Writing...\n";
	out->Write(ref_seq.get(), size);

	return true;
}

// ************************************************************************************
// Compact direct reference sequence (2 symbols -> 1 byte)
bool ref_seq_dir_compact(const string index_name, const string dest_dir, const string temp_dir)
{
	cout << "\n***** Stage 3. Reference sequence compaction (direct)\n";

	shared_ptr<CMapperFile> in(new CMapperFile(EXT_REF_SEQ_DIR, MARKER_REF_SEQ_DIR));
	if(!in->OpenRead(temp_dir + index_name))
	{
		cout << "Cannot open reference sequence (direct)\n";
		return false;
	}
	int64_t size = in->GetSize();

	shared_ptr<CMapperFile> out(new CMapperFile(EXT_REF_SEQ_DIR_PCK, MARKER_REF_SEQ_DIR_PCK));
	if(!out->OpenWrite(dest_dir + index_name))
	{
		cout << "Cannot create compacted reference sequence (direct)\n";
		return false;
	}

	cout << "Reading...\n";
	shared_ptr<char> ref_seq(new char[size+1]);
	in->Read(ref_seq.get(), size);
	ref_seq.get()[size] = sym_code_N_ref;

	cout << "Compacting...\n";
	char *ptr = ref_seq.get();
	for(int64_t i = 0; i < (size+1)/2; ++i)
		ptr[i] = (ptr[2*i] << 4) + ptr[2*i+1];

	cout << "Writing...\n";
	out->Write(ref_seq.get(), (size+1)/2);
	
	return true;
}

// ************************************************************************************
// Compact rev. comp. reference sequence (2 symbols -> 1 byte)
bool ref_seq_rc_compact(const string index_name, const string dest_dir, const string temp_dir)
{
	cout << "\n***** Stage 4. Reference sequence compaction (rev. comp.)\n";

	shared_ptr<CMapperFile> in(new CMapperFile(EXT_REF_SEQ_RC, MARKER_REF_SEQ_RC));
	if(!in->OpenRead(temp_dir + index_name))
	{
		cout << "Cannot open reference sequence (rev. comp.)\n";
		return false;
	}
	int64_t size = in->GetSize();

	shared_ptr<CMapperFile> out(new CMapperFile(EXT_REF_SEQ_RC_PCK, MARKER_REF_SEQ_RC_PCK));
	if(!out->OpenWrite(dest_dir + index_name))
	{
		cout << "Cannot create compacted reference sequence (rev. comp.)\n";
		return false;
	}

	cout << "Reading...\n";
	shared_ptr<char> ref_seq(new char[size+1]);
	in->Read(ref_seq.get(), size);
	ref_seq.get()[size] = sym_code_N_ref;

	cout << "Compacting...\n";
	char *ptr = ref_seq.get();
	for(int64_t i = 0; i < (size+1)/2; ++i)
		ptr[i] = (ptr[2*i] << 4) + ptr[2*i+1];

	cout << "Writing...\n";
	out->Write(ref_seq.get(), (size+1)/2);

	return true;
}

// EOF
