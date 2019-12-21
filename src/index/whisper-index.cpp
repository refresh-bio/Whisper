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


#include <cstdio>
#include "../common/defs.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "ref_seq.h"
#include "suf_arr.h"
#include "vcf.h"

using namespace std;

string index_name;
vector<string> ref_seq_names;
vector<string> vcf_names;
string dest_dir;
string temp_dir;
int stage_no;

void usage();
bool parse_params(int argc, char **argv);
bool remove_temps(string index_name, string temp_dir);

// ************************************************************************************
void usage()
{
	cerr << "Whisper index construction v. " << MAPPER_VERSION << "\n";
	cerr << "Usage: \n";
#ifdef _DEVELOPMENT_MODE
	cerr << "   whisper-index <index_name> <ref_seq_file_name> <dest_dir> <temp_dir> <vcf_name> [stage_no]\n";
	cerr << "   whisper-index <index_name> <@ref_seq_files_name> <dest_dir> <temp_dir> <@vcf_files_name> [stage_no]\n";
	cerr << "Stages:\n";
	cerr << "   0 - all stages (default)\n";
	cerr << "   1 - ref. seq. construction (direct)\n";
	cerr << "   2 - ref. seq. construction (rev. comp.)\n";
	cerr << "   3 - ref. seq. compaction (direct)\n";
	cerr << "   4 - ref. seq. compaction (rev. comp.)\n";
	cerr << "   5 - suffix array construction (direct)\n";
	cerr << "   6 - suffix array construction (rev. comp.)\n";
	cerr << "   7 - LUT for suffix array construction (direct)\n";
	cerr << "   8 - LUT for suffix array construction (rev. comp.)\n";
	cerr << "   9 - SNP and short indel construction from VCF file\n";
	cerr << "   99 - removing all temporary files\n";
#else
	cerr << "   whisper-index <index_name> <ref_seq_file_name> <dest_dir> <temp_dir>\n";
	cerr << "   whisper-index <index_name> <@ref_seq_files_name> <dest_dir> <temp_dir>\n";
#endif
	cerr << "Hints:\n";
	cerr << "   * vcf_name can be . if not used\n";

	exit(0);
}

// ************************************************************************************
bool parse_params(int argc, char **argv)
{
	index_name = string(argv[1]);
	dest_dir   = NormalizeDirectory(string(argv[3]));
	temp_dir   = NormalizeDirectory(string(argv[4]));
//	dbsnp_name = string(argv[5]);

	if(argv[2][0] != '@')
		ref_seq_names.push_back(string(argv[2]));
	else
	{
		ifstream inf(argv[2]+1);
		if(!inf.good())
		{
			cerr << "Cannot open file : " << string(argv[2]+1) << "\n";
			return false;
		}

		string s;
		while(getline(inf, s))
			if(s != "")
				ref_seq_names.push_back(s);
		inf.close();
	}

	if (argv[5][0] != '@')
		vcf_names.push_back(string(argv[5]));
	else
	{
		ifstream inf(argv[5] + 1);
		if (!inf.good())
		{
			cerr << "Cannot open file : " << string(argv[5] + 1) << "\n";
			return false;
		}

		string s;
		while (getline(inf, s))
			if (s != "")
				vcf_names.push_back(s);
		inf.close();
	}

#ifdef _DEVELOPMENT_MODE
	if(argc == 7)
		stage_no = atoi(argv[6]);
	else
#endif
		stage_no = 0;
		
	return !ref_seq_names.empty();
}

// ************************************************************************************
bool remove_temps(string index_name, string temp_dir)
{
	cerr << "\n***** Stage 9. Removing temporary files\n";

	string name = temp_dir + index_name + EXT_REF_SEQ_DIR;
	cerr << "Removing " << name << "\n";
	remove(name.c_str());

	name = temp_dir + index_name + EXT_REF_SEQ_RC;
	cerr << "Removing " << name << "\n";
	remove(name.c_str());

	return true;
}

// ************************************************************************************
int main(int argc, char **argv)
{
	if(argc < 5 || !parse_params(argc, argv))
		usage();

	// Reference sequence (direct) construction
	if(stage_no == 0 || stage_no == 1)
		if(!ref_seq_dir_construct(index_name, ref_seq_names, dest_dir, temp_dir))
			return 0;

	// Reference sequence (rev. comp.) construction
	if(stage_no == 0 || stage_no == 2)
		if(!ref_seq_rc_construct(index_name, temp_dir))
			return 0;

	// Reference sequence (direct) compaction
	if(stage_no == 0 || stage_no == 3)
		if(!ref_seq_dir_compact(index_name, dest_dir, temp_dir))
			return 0;

	// Reference sequence (rev. comp.) compation
	if(stage_no == 0 || stage_no == 4)
		if(!ref_seq_rc_compact(index_name, dest_dir, temp_dir))
			return 0;

	// Suffix array (direct) construction
	if(stage_no == 0 || stage_no == 5)
		if(!suffix_array_construct_dir(index_name, dest_dir, temp_dir))
			return 0;

	// Suffix array (rev. comp.) construction
	if(stage_no == 0 || stage_no == 6)
		if(!suffix_array_construct_rc(index_name, dest_dir, temp_dir))
			return 0;

	// LUT (direct) computation
	if(stage_no == 0 || stage_no == 7)
		if(!lut_compute_dir(index_name, dest_dir, temp_dir))
			return 0;

	// LUT (rev. comp.) computation
	if(stage_no == 0 || stage_no == 8)
		if(!lut_compute_rc(index_name, dest_dir, temp_dir))
			return 0;

	// SNP and short variants construction
#ifdef ENABLE_VCF_VARIANTS
	if (stage_no == 0 || stage_no == 9)
		if (!vcf_construct(index_name, dest_dir, vcf_names))
			return 0;
#endif

	// Remove temporary files
	if(stage_no == 0 || stage_no == 99)
		if(!remove_temps(index_name, temp_dir))
			return 0;
			
	cerr << "\nCompleted.\n";
	
	return 0;
}

// EOF
