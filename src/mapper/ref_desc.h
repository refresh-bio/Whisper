// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.0
// Date    : 2017-11-30
// License : GNU GPL 3
// *******************************************************************************************


#ifndef _REF_DESC_H
#define _REF_DESC_H

#include "../common/defs.h"
#include "../common/idstore.h"
#include "../common/queue.h"
#include "../common/mmgr.h"
#include "../common/stats.h"
#include "../common/types.h"
#include "../common/params.h"
#include "../common/joiner_mgr.h"

#include <vector>

using namespace std;

// ************************************************************************************
class CRefSeqDesc
{
	vector<seq_desc_t> seq_desc;
	string file_name;

public:
	CRefSeqDesc();
	~CRefSeqDesc();

	bool Load(string _file_name);

	bool Translate(uint32_t raw_pos, string &seq_name, int32_t &pos, int32_t& id, const uchar_t* msg);
	bool GetDescription(vector<seq_desc_t> &_seq_desc);
	vector<seq_desc_t>& GetDescription() { return seq_desc;  }

	const seq_desc_t& operator[](int i) const { return seq_desc[i]; }
	seq_desc_t& operator[](int i) { return seq_desc[i]; }

	int64_t GetSize(string _file_name);
};

#endif

// EOF
