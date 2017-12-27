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


#include "sa.h"
#include <iostream>
#include <thread>

// ************************************************************************************
//
CSuffixArray::CSuffixArray(bool _cache_sa)
{
	is_valid = false;

	file_sa		   = nullptr;
	file_lut_short = nullptr;
	file_lut_long  = nullptr;

	data_lut       = nullptr;

	cache_sa       = _cache_sa;
}

// ************************************************************************************
//
CSuffixArray::~CSuffixArray()
{
	close_files();

	for(auto &p: parts)
		delete[] p.second.first;
}

// ************************************************************************************
// 1. Set input file names
// 2. Open the files
bool CSuffixArray::SetIndexName(string _index_name, genome_t genome_type)
{
	unique_lock<mutex> lck(mtx);

	close_files();
	if(cache_sa)
		for(auto &p: parts)
			delete[] p.second.first;

	string sa_ext, lut_short_ext, lut_long_ext;
	string sa_marker, lut_short_marker, lut_long_marker;

	if(genome_type == genome_t::direct)
	{
		sa_ext           = EXT_SA_DIR;
		lut_short_ext    = EXT_LUT_SHORT_DIR;
		lut_long_ext     = EXT_LUT_LONG_DIR;
		sa_marker        = MARKER_SA_DIR;
		lut_short_marker = MARKER_LUT_SHORT_DIR;
		lut_long_marker  = MARKER_LUT_LONG_DIR;
	}
	else
	{
		sa_ext           = EXT_SA_RC;
		lut_short_ext    = EXT_LUT_SHORT_RC;
		lut_long_ext     = EXT_LUT_LONG_RC;
		sa_marker        = MARKER_SA_RC;
		lut_short_marker = MARKER_LUT_SHORT_RC;
		lut_long_marker  = MARKER_LUT_LONG_RC;
	}

	string dir_name = _index_name;

	file_sa.reset(new CMapperFile(sa_ext, sa_marker, sizeof(uint32_t)));
	if(!file_sa->OpenRead(dir_name))
		return false;

	file_lut_short.reset(new CMapperFile(lut_short_ext, lut_short_marker, sizeof(uint32_t)));
	if(!file_lut_short->OpenRead(dir_name))
	{
		file_sa.reset();
		return false;
	}
	size_lut_short = file_lut_short->GetSize();

	file_lut_long.reset(new CMapperFile(lut_long_ext, lut_long_marker, sizeof(uint32_t)));
	if(!file_lut_long)
	{
		file_sa.reset();
		file_lut_short.reset();
		return false;
	}

	data_lut = new uint32_t[size_lut_short];
	if(file_lut_short->Read(data_lut, size_lut_short) != size_lut_short)
	{
		close_files();
		return false;
	}
	file_lut_short.reset();

	is_valid = true;

	file_sa->Seek(0);
	file_sa->Read(&ref_seq_size, 1);

	return true;
}

// ************************************************************************************
// Return size of the reference sequence
uint32_t CSuffixArray::GetRefSize()
{
	if(!is_valid)
		return 0;

	return ref_seq_size;
}

// ************************************************************************************
//
void CSuffixArray::close_files()
{
	if(!is_valid)
		return;

	file_sa.reset();
	file_lut_short.reset();
	file_lut_long.reset();

	if(data_lut)
	{
		delete[] data_lut;
		data_lut = nullptr;
	}

	is_valid = false;
}

// ************************************************************************************
// Return the part of suffix array
// The client must release the memory for data
bool CSuffixArray::GetSAPart(uint32_t prefix, uint32_t prefix_len, uint32_t* &data, uint32_t &size)
{
	unique_lock<mutex> lck(mtx);
	
	if(!is_valid)
		return false;

	if(prefix_len > lut_short_prefix_len)
		return false;

	uint32_t filled_prefix      = prefix     << (2 * (lut_short_prefix_len - prefix_len));
	uint32_t filled_prefix_next = (prefix+1) << (2 * (lut_short_prefix_len - prefix_len));

	uint32_t first_index = data_lut[filled_prefix];
	uint32_t last_index  = data_lut[filled_prefix_next];

	size = last_index - first_index;

	if(cache_sa)
	{
		pair<uint32_t, uint32_t> idx = make_pair(prefix, prefix_len);
		auto p = parts.find(idx);
		if(p == parts.end())
		{
			data = new uint32_t[size];
			file_sa->Seek(first_index);
			file_sa->Read(data, size);
			parts[idx] = make_pair(data, size);
		}
		else
		{
			data = p->second.first;
			size = p->second.second;
		}
	}
	else
	{
		data = new uint32_t[size];
		file_sa->Seek(first_index);
		file_sa->Read(data, size);
	}

	return true;
}

// ************************************************************************************
bool CSuffixArray::ReleaseSAPart(uint32_t* &data)
{
//	unique_lock<mutex> lck(mtx);
	
	if(!is_valid)
		return false;

	if(!cache_sa && data)
		delete[] data;

	return true;
}

// ************************************************************************************
// Return size of the part of suffix array of given prefix
uint32_t CSuffixArray::GetSAPartSize(uint32_t prefix, uint32_t prefix_len)
{
	if(!is_valid)
		return false;

	if(prefix_len > lut_short_prefix_len)
		return false;

	uint32_t filled_prefix      = prefix     << (2 * (lut_short_prefix_len - prefix_len));
	uint32_t filled_prefix_next = (prefix+1) << (2 * (lut_short_prefix_len - prefix_len));

	uint32_t first_index = data_lut[filled_prefix];
	uint32_t last_index  = data_lut[filled_prefix_next];

	uint32_t size = last_index - first_index;

	return size;
}

// ************************************************************************************
// Return LUT (short)
bool CSuffixArray::GetShortLUT(uint32_t* &_data_lut, uint32_t &_lut_short_prefix_len)
{
	if(!is_valid)
		return false;

	_data_lut			  = data_lut;
	_lut_short_prefix_len = lut_short_prefix_len;

	return true;
}

// EOF
