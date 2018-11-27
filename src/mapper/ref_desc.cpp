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


#include "ref_desc.h"
#include <algorithm>

// ************************************************************************************
CRefSeqDesc::CRefSeqDesc()
{
}

// ************************************************************************************
CRefSeqDesc::~CRefSeqDesc()
{
}

// ************************************************************************************
bool CRefSeqDesc::Load(string _file_name)
{
	uint64_t seq_size, pos_in_ref_seq, no_initial_Ns, no_final_Ns, no_starting_Ns;
	uint32_t seq_name_len;

	file_name = _file_name;

	shared_ptr<CMapperFile> in_desc(new CMapperFile(EXT_REF_SEQ_DESC, MARKER_REF_SEQ_DESC, 1));
	if(!in_desc->OpenRead(file_name))
		return false;

	seq_desc.clear();

	int32_t id = 0;
	while(!in_desc->Eof())
	{
		in_desc->Read(&seq_name_len, sizeof(uint32_t));		// seq. name length
		string seq_name;
		bool trimming = false;
		for(uint32_t i = 0; i < seq_name_len; ++i)
		{
			uchar_t c;
			in_desc->Read(&c, 1);

			if (!trimming)
			{
				if (c == ' ' || c == '\t' || c == '\r' || c == '\n')
					trimming = true;

				if (!trimming)
					seq_name.push_back(c);
			}
		}
		in_desc->Read(&seq_size, sizeof(uint64_t));
		in_desc->Read(&pos_in_ref_seq, sizeof(uint64_t));
		in_desc->Read(&no_initial_Ns, sizeof(uint64_t));
		in_desc->Read(&no_final_Ns, sizeof(uint64_t));
		in_desc->Read(&no_starting_Ns, sizeof(uint64_t));

		seq_desc.push_back(seq_desc_t(seq_name, id++, seq_size, pos_in_ref_seq, no_initial_Ns, no_final_Ns, no_starting_Ns));
	}

	return true;
}

// ************************************************************************************
// Return reference size
int64_t CRefSeqDesc::GetSize(string _file_name)
{
	return FileSize(_file_name + EXT_REF_SEQ_DESC);
}

// ************************************************************************************
// Translate raw seq. position into pair: (reference seq. name, position in this sequence)
bool CRefSeqDesc::Translate(uint32_t raw_pos, string &seq_name, int32_t &pos, int32_t& id, const uchar_t* msg)
{
	uint32_t raw_pos1 = raw_pos + 1;
	vector<seq_desc_t>::iterator p = lower_bound(seq_desc.begin(), seq_desc.end(), raw_pos1, 
		[](const seq_desc_t& x, uint32_t val)->bool {return x.pos_in_ref_seq + SEPARATE_N_LEN / 2 < val; });	// to support clipped positions that can be prior to first valid base in reference

	if(p == seq_desc.end())
	{
		if(raw_pos >= seq_desc.back().pos_in_ref_seq + seq_desc.back().size + 2 * SEPARATE_N_LEN)
		{
			seq_name = "";
			pos	= 0;
	
			cerr << "Cannot find raw_pos: " << raw_pos << ", info: " << msg << endl;
			return false;
		}
	}

	if(p == seq_desc.begin())
	{
		cerr << "Error in Translate! raw_pos1: " << raw_pos1 << ", info:" << msg << endl;
		return false;
	}
	--p;

	seq_name = p->name;
	id = p->id;
	
	int64_t r = raw_pos;
	r -= p->pos_in_ref_seq;
	r += p->no_initial_Ns;
	r -= p->no_starting_Ns;
	++r;						// convert from 0-base to 1-base

	pos = (int32_t) r;

	return true;
}

// ************************************************************************************
// Return description of sequences in the collection
bool CRefSeqDesc::GetDescription(vector<seq_desc_t> &_seq_desc)
{
	_seq_desc = seq_desc;

	return !seq_desc.empty();
}

// EOF
