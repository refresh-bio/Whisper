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


#ifndef _IDSTORE_H
#define _IDSTORE_H

#include "defs.h"
#include <string>
#include <map>
#include <utility>
#include <memory>
#include <mutex>
#include <cstdio>

#include "../common/utils.h"

using namespace std;


// ************************************************************************************
// Map blocks of FASTQ files to unique ids
// ************************************************************************************
class CIDStore
{
	mutable mutex mtx;

	map<pair<uint32_t, uint32_t>, pair<string, read_id_t>> mappings;
	uint32_t total_bits;
	uint32_t sub_block_bits;
	uint32_t in_block_bits;

	uint32_t mappings_effective_size;
	uint32_t verbosity_level;

public:
	CIDStore(uint32_t _total_bits, uint32_t _sub_block_bits, uint32_t _in_block_bits, uint32_t _verbosity_level);
	~CIDStore();

	bool Save(shared_ptr<CMapperFile> file);
	bool Load(shared_ptr<CMapperFile> file);

	read_id_t RegisterBlock(uint32_t file_no, string file_name, uint32_t file_part, read_id_t prev_id);
	read_id_t GetBlockID(uint32_t file_no, string file_name, uint32_t file_part);
	uint32_t GetNumGroups();

	bool GetIDBits(uint32_t &_total_bits, uint32_t &_sub_block_bits, uint32_t &_in_block_bits);
};

#endif

// EOF
