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


#ifndef _REFERENCE_H
#define _REFERENCE_H

#include <string>
#include <cstdio>
#include <memory>
#include <mutex>
#include "../common/defs.h"
#include "../common/utils.h"

using namespace std;

// ************************************************************************************
class CReference
{
	uint32_t size;
	uchar_t *data;

	mutable mutex mtx;

public:
	CReference();
	~CReference();

	bool SetIndexName(string _index_name);

	uchar_t *GetData();
	uint32_t GetSize();
};

#endif

// EOF
