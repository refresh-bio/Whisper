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


#include "reference.h"

// ************************************************************************************
//
CReference::CReference()
{
	data = nullptr;
}

// ************************************************************************************
CReference::~CReference()
{
	if(data)
		delete[] data;
}

// ************************************************************************************
// 1. Set file name of the reference sequence
// 2. Open this file
// 3. Read the data
bool CReference::SetIndexName(string _index_name)
{
	unique_lock<mutex> lck(mtx);

	shared_ptr<CMapperFile> file(new CMapperFile(EXT_REF_SEQ_DIR_PCK, MARKER_REF_SEQ_DIR_PCK));

	if(data)
	{
		delete[] data;
		data = nullptr;
	}

	// Try to open and check the file format
	if(!file->OpenRead(_index_name))
		return false;
	size = (uint32_t) file->GetSize();

	// Read the data
	data = new uchar_t[size];
	file->Read(data, size);
	
	return true;
}

// ************************************************************************************
//
uchar_t *CReference::GetData()
{
	unique_lock<mutex> lck(mtx);

	return data;
}

// ************************************************************************************
//
uint32_t CReference::GetSize()
{
	unique_lock<mutex> lck(mtx);

	return size;
}

// EOF
