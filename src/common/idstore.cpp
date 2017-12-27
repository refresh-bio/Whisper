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


#include "idstore.h"
#include <iostream>
#include <memory>

// ************************************************************************************
// CIDStore
// ************************************************************************************

// ************************************************************************************
CIDStore::CIDStore(uint32_t _total_bits, uint32_t _sub_block_bits, uint32_t _in_block_bits, uint32_t _verbosity_level)
{
	total_bits     = _total_bits;
	sub_block_bits = _sub_block_bits;
	in_block_bits  = _in_block_bits;

	mappings_effective_size = 0;

	verbosity_level = _verbosity_level;
}

// ************************************************************************************
CIDStore::~CIDStore()
{
}

// ************************************************************************************
// Save mappings data (file_name, part_in_file) - id to file
bool CIDStore::Save(shared_ptr<CMapperFile> file)
{
	unique_lock<mutex> lck(mtx);

	if(mappings.empty())
		return false;

	file->Write(&total_bits, sizeof(uint32_t));
	file->Write(&sub_block_bits, sizeof(uint32_t));
	file->Write(&in_block_bits, sizeof(uint32_t));
	file->Write(&mappings_effective_size, sizeof(uint32_t));

	for(auto &p : mappings)
	{
		file->Write(p.second.first.c_str(), p.second.first.length()+1);	// file_name
		file->Write(&(p.first.first), sizeof(uint32_t));				// file_no
		file->Write(&(p.first.second), sizeof(uint32_t));				// part_no
		file->Write(&(p.second.second), sizeof(read_id_t));				// block_id
	}

	return true;
}

// ************************************************************************************
// Load mappings data from file
bool CIDStore::Load(shared_ptr<CMapperFile> file)
{
	unique_lock<mutex> lck(mtx);

	int64_t size = file->GetSize();
	string data;

	if(!size)
		return false;

	file->Read(&total_bits, sizeof(uint32_t));
	file->Read(&sub_block_bits, sizeof(uint32_t));
	file->Read(&in_block_bits, sizeof(uint32_t));
	file->Read(&mappings_effective_size, sizeof(uint32_t));

	mappings.clear();

	while(file->ReadString(data))
	{
		uint32_t file_no, part_no;
		read_id_t block_id;
		file->Read(&file_no, sizeof(uint32_t));
		file->Read(&part_no, sizeof(uint32_t));
		file->Read(&block_id, sizeof(read_id_t));
		mappings[make_pair(file_no, part_no)] = make_pair(data, block_id);
	}

	return true;
}

// ************************************************************************************
// Return unique ID range for block of FASTQ files
// The range is of size (1 << in_block_bits)
read_id_t CIDStore::RegisterBlock(uint32_t file_no, string file_name, uint32_t file_part, read_id_t prev_id)
{
	unique_lock<mutex> lck(mtx);
	read_id_t id = prev_id;
	read_id_t sub_block_mask = (((read_id_t) 1) << sub_block_bits) - 1;

	if(((prev_id >> in_block_bits) & sub_block_mask) == sub_block_mask)	// new block_id necessary
	{
		auto p = mappings.find(make_pair(file_no ^ 1, file_part));			// look for paired FASTQ block
		if(p == mappings.end())
			id = (((read_id_t) (mappings_effective_size++)) << (in_block_bits + sub_block_bits)) + (file_no & 1);
		else
			id = p->second.second ^ 1;							// block_id of the paired block
		mappings[make_pair(file_no, file_part)] = make_pair(file_name, id);
	}
	else
	{
		id  = prev_id;
		id += ((read_id_t) 1) << in_block_bits;
		id &= ~(((read_id_t) 1 << (in_block_bits)) - 1);
		id += file_no & 1;
	}

	if(verbosity_level > 2)
		cout << "IDStore:-register " << hex << id << "\n";

	return id;
}

// ************************************************************************************
// Return id of block (when all ids are already assigned)
read_id_t CIDStore::GetBlockID(uint32_t file_no, string file_name, uint32_t file_part)
{
	unique_lock<mutex> lck(mtx);
	read_id_t id;
	read_id_t sub_block_mask = (((read_id_t) 1) << sub_block_bits) - 1;

	auto p = mappings.find(make_pair(file_no, file_part & ~sub_block_mask));
	if(p == mappings.end())
		return ~((read_id_t) 0);

	id = p->second.second;
	id += ((read_id_t) file_part & sub_block_mask) << in_block_bits;

	if(verbosity_level > 2)
		cout << "IDStore-get: " << hex << id << "\n";

	return id;
}

// ************************************************************************************
// Return number of groups
uint32_t CIDStore::GetNumGroups()
{
	unique_lock<mutex> lck(mtx);

	return mappings_effective_size;
}

// ************************************************************************************
bool CIDStore::GetIDBits(uint32_t &_total_bits, uint32_t &_sub_block_bits, uint32_t &_in_block_bits)
{
	unique_lock<mutex> lck(mtx);

	_total_bits     = total_bits;
	_sub_block_bits = sub_block_bits;
	_in_block_bits  = in_block_bits;

	return true;
}

// EOF
