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


#ifndef _REF_SEQ_H
#define _REF_SEQ_H

#include "../common/defs.h"
#include "../common/types.h"
#include <string>
#include <vector>

using namespace std;

bool ref_seq_append(char *seq, uint64_t &pos, const string file_name);

bool ref_seq_dir_construct(const string index_name, const vector<string> &ref_seq_names, const string dest_dir, const string temp_dir);
bool ref_seq_rc_construct(const string index_name, const string temp_dir);

bool ref_seq_dir_compact(const string index_name, const string dest_dir, const string temp_dir);
bool ref_seq_rc_compact(const string index_name, const string dest_dir, const string temp_dir);

#endif

// EOF
