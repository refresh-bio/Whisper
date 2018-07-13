// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.1
// Date    : 2018-07-10
// License : GNU GPL 3
// *******************************************************************************************


#ifndef _SUF_ARR_H
#define _SUF_ARR_H

#include <string>
#include <memory>
#include "../common/defs.h"
#include "../common/utils.h"

using namespace std;

bool suffix_array_construct_dir(const string index_name, const string dest_dir, const string temp_dir);
bool suffix_array_construct_rc(const string index_name, const string dest_dir, const string temp_dir);

bool lut_compute_dir(const string index_name, const string dest_dir, const string temp_dir);
bool lut_compute_rc(const string index_name, const string dest_dir, const string temp_dir);

bool suffix_array_construct(shared_ptr<CMapperFile> in, shared_ptr<CMapperFile> out, uint32_t size, bool translate_Ns);
bool suffix_array_cleanup(shared_ptr<uchar_t> text, shared_ptr<uint32_t> sa, uint32_t &size, uint32_t prefix_len);

bool lut_compute(shared_ptr<CMapperFile> in_text, shared_ptr<CMapperFile> in_sa, shared_ptr<CMapperFile> out_short, shared_ptr<CMapperFile> out_long, uint32_t size_text, uint32_t size_sa);

#endif

// EOF
