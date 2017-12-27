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


#ifndef _SA_H
#define _SA_H

#include <string>
#include <cstdio>
#include <memory>
#include <mutex>
#include <map>
#include <utility>

#include "../common/defs.h"
#include "../common/utils.h"

using namespace std;

// ************************************************************************************
class CSuffixArray
{
	bool is_valid;

	mutable mutex mtx;

	shared_ptr<CMapperFile> file_sa;
	shared_ptr<CMapperFile> file_lut_short;
	shared_ptr<CMapperFile> file_lut_long;

	map<pair<uint32_t, uint32_t>, pair<uint32_t*, uint32_t>> parts;

	uint32_t *data_lut;
	uint32_t ref_seq_size;
	int64_t size_lut_short;
	bool cache_sa;

	void close_files();

public:
	CSuffixArray(bool _cache_sa = false);
	~CSuffixArray();

	bool SetIndexName(string _index_name, genome_t genome_type);
	uint32_t GetRefSize();

	bool GetSAPart(uint32_t prefix, uint32_t prefix_len, uint32_t* &data, uint32_t &size);
	bool ReleaseSAPart(uint32_t* &data);

	bool GetShortLUT(uint32_t* &_data_lut, uint32_t &_lut_short_prefix_len);

	uint32_t GetSAPartSize(uint32_t prefix, uint32_t prefix_len);
};

#endif

// EOF
