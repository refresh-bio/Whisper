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


#include "utils.h"
#include <malloc.h>
#include <iostream>
#include <memory>
#include <sstream>

// ************************************************************************************
void *ptr_align(void *ptr, size_t alignment)
{
	// Check whether alignment is a power of 2
	if (alignment & (alignment - 1))
		return nullptr;

	uint64_t iptr = reinterpret_cast<uint64_t>(ptr);
	uint64_t shift = alignment - iptr % alignment;

	// Shift pointer as necessary
	return reinterpret_cast<void*>(reinterpret_cast<char*>(ptr) + shift);
}

// ************************************************************************************
void *alloc_aligned(void *&raw_ptr, size_t size, size_t alignment)
{
	raw_ptr = realloc(raw_ptr, size + alignment - 1);

	uint64_t iptr = reinterpret_cast<uint64_t>(raw_ptr);
	uint64_t shift = alignment - iptr % alignment;

	// Shift pointer as necessary
	return reinterpret_cast<void*>(reinterpret_cast<char*>(raw_ptr) + shift);
}


// ************************************************************************************
// CMapperFile
// ************************************************************************************

// ************************************************************************************
CMapperFile::CMapperFile()
{
	state = open_mode_t::not_set;
	SetParams("", "", 1);
}
	
// ************************************************************************************
CMapperFile::CMapperFile(string _extension, string _marker, uint32_t _element_size)
{
	state = open_mode_t::not_set;
	SetParams(_extension, _marker, _element_size);
}
	
// ************************************************************************************
CMapperFile::~CMapperFile()
{
	if(state == open_mode_t::opened_for_read || state == open_mode_t::opened_for_write)
		Close();
}

// ************************************************************************************
bool CMapperFile::SetParams(string _extension, string _marker, uint32_t _element_size)
{
	if(state != open_mode_t::closed && state != open_mode_t::not_set)
		return false;

	name         = "";
	extension    = _extension;
	marker       = _marker;
	element_size = _element_size;

	if(element_size == 0)
		element_size = 1;
	marker_length = marker.length();

	state        = open_mode_t::closed;

	return true;
}

// ************************************************************************************
bool CMapperFile::OpenRead(string _name)
{
	if(state != open_mode_t::closed && state != open_mode_t::not_set)
		return false;

	name = _name;
	if(extension != "")
		name += extension;

	file = fopen(name.c_str(), "rb");

	if(!file)
		return file_error("File " + name + " does not exist\n");

	int32_t marker_length = (int32_t) marker.length();
	shared_ptr<char> tmp(new char[marker_length+1]);

	if(marker_length)
	{
		if((int32_t) fread(tmp.get(), 1, marker_length, file) != marker_length)
			return file_error("Wrong marker in " + name + " file\n");
		if(memcmp(tmp.get(), marker.c_str(), marker_length) != 0)
			return file_error("Wrong marker in " + name + " file\n");
	}

	seek(-marker_length, SEEK_END);

	size = (tell() - marker_length) / element_size;

	if(marker_length)
	{
		if((int32_t) fread(tmp.get(), 1, marker_length, file) != marker_length)
			return file_error("Wrong marker in " + name + " file\n");
		if(memcmp(tmp.get(), marker.c_str(), marker_length) != 0)
			return file_error("Wrong marker in " + name + " file\n");
	}

	seek(marker_length, SEEK_SET);

	state = open_mode_t::opened_for_read;

	return true;
}
	
// ************************************************************************************
bool CMapperFile::OpenWrite(string _name)
{
	if(state != open_mode_t::closed && state != open_mode_t::not_set)
		return false;

	name = _name;
	if(extension != "")
		name += extension;

	file = fopen(name.c_str(), "wb");

	if(!file)
		file_error("Cannot create file: " + name + "\n");

	if(marker_length)
		fwrite(marker.c_str(), 1, marker_length, file);
	
	state = open_mode_t::opened_for_write;

	return true;
}
	
// ************************************************************************************
bool CMapperFile::OpenStream(FILE *stream)
{
	file = stream;
	state = open_mode_t::direct_stream;

	return true;
}

// ************************************************************************************
bool CMapperFile::Close()
{
	if(state == open_mode_t::opened_for_read)
	{
		fclose(file);
		state = open_mode_t::closed;
		return true;
	}

	if(state == open_mode_t::opened_for_write)
	{
		if(marker_length)
			fwrite(marker.c_str(), 1, marker_length, file);
		fclose(file);
		state = open_mode_t::closed;
		return true;
	}

	return false;
}

// ************************************************************************************
bool CMapperFile::Remove()
{
	if(state != open_mode_t::closed)
		return false;

	remove(name.c_str());

	state = open_mode_t::not_set;

	return true;
}

// ************************************************************************************
int64_t CMapperFile::GetSize()
{
	if(state == open_mode_t::opened_for_read)
		return size;

	return 0;
}
	
// ************************************************************************************
int64_t CMapperFile::GetRawSize()
{
	if(state == open_mode_t::opened_for_read)
		return size * element_size;

	return 0;
}
	
// ************************************************************************************
int64_t CMapperFile::Read(void* data, int64_t count)
{
	if(state != open_mode_t::opened_for_read)
		return -1;

	return fread(data, element_size, count, file);
}

// ************************************************************************************
int64_t CMapperFile::ReadUint16(uint16_t &data)
{
	if (state != open_mode_t::opened_for_read && state != open_mode_t::direct_stream)
		return -1;

	data = 0;

	for (int i = 0; i < 2; ++i)
	{
		int x = getc(file);

		if (x == EOF)
			return -1;

		data += ((uint16_t)x) << (8 * i);
	}

	return 2;
}

// ************************************************************************************
int64_t CMapperFile::ReadUint32(uint32_t &data)
{
	if (state != open_mode_t::opened_for_read && state != open_mode_t::direct_stream)
		return -1;

	data = 0;

	for (int i = 0; i < 4; ++i)
	{
		int x = getc(file);

		if (x == EOF)
			return -1;

		data += ((uint32_t)x) << (8 * i);
	}

	return 4;
}

// ************************************************************************************
int64_t CMapperFile::ReadString(string &data)
{
	data = "";
	while(!Eof())
	{
		int c = getc(file);
		if(c)
			data.push_back((char) c);
		else
			break;
	}

	return data.length();
}

// ************************************************************************************
int64_t CMapperFile::Write(void* data, int64_t count)
{
	if(state != open_mode_t::opened_for_write && state != open_mode_t::direct_stream)
		return -1;

	return fwrite(data, element_size, count, file);
}
	
// ************************************************************************************
int64_t CMapperFile::Write(const void* data, int64_t count)
{
	if(state != open_mode_t::opened_for_write && state != open_mode_t::direct_stream)
		return -1;

	return fwrite(data, element_size, count, file);
}
	
// ************************************************************************************
int64_t CMapperFile::WriteByte(uint32_t data)
{
	if (state != open_mode_t::opened_for_write && state != open_mode_t::direct_stream)
		return -1;

	putc(data, file);

	return 1;
}

// ************************************************************************************
int64_t CMapperFile::WriteUint16(uint16_t data)
{
	if (state != open_mode_t::opened_for_write && state != open_mode_t::direct_stream)
		return -1;

	for (int i = 0; i < 2; ++i)
	{
		auto x = data & 0xff;
		data >>= 8;

		putc(x, file);
	}

	return 2;
}

// ************************************************************************************
int64_t CMapperFile::WriteUint32(uint32_t data)
{
	if (state != open_mode_t::opened_for_write && state != open_mode_t::direct_stream)
		return -1;

	for (int i = 0; i < 4; ++i)
	{
		auto x = data & 0xff;
		data >>= 8;

		putc(x, file);
	}

	return 4;
}

// ************************************************************************************
int64_t CMapperFile::SeekBegin()
{
	if(state != open_mode_t::opened_for_read)
		return -1;

	return seek(marker_length, SEEK_SET);
}
	
// ************************************************************************************
int64_t CMapperFile::Seek(int64_t pos)
{
	if(state != open_mode_t::opened_for_read)
		return -1;

	return seek(marker_length + pos * element_size, SEEK_SET);
}
	
// ************************************************************************************
bool CMapperFile::Eof()
{
	if(!file)
		return false;

	return feof(file) || tell() >= size * element_size + (int64_t) marker_length;
}

// ************************************************************************************
void CMapperFile::SetBufSize(uint64_t size)
{
	if(state != open_mode_t::closed && state != open_mode_t::direct_stream)
		setvbuf(file, nullptr, _IOFBF, size);
}

// ************************************************************************************
bool CMapperFile::file_error(string error_code)
{
	cerr << error_code;
	return false;
}

// ************************************************************************************
inline int64_t CMapperFile::seek(int64_t offset, int origin)
{
#ifdef WIN32
	return _fseeki64(file, offset, origin);
#else
	return fseeko64(file, offset, origin);
#endif
}

// ************************************************************************************
inline int64_t CMapperFile::tell()
{
#ifdef WIN32
	return _ftelli64(file);
#else
	return ftello64(file);
#endif
}


// ************************************************************************************
// CNumericConversions statics
uchar_t CNumericConversions::digits[];
CNumericConversions::_si CNumericConversions::_init;
uint64_t CNumericConversions::powers10[];


// ************************************************************************************
void StoreUIntLSB(vector<uchar_t> &dest, uint64_t data, uint32_t no_bytes)
{
	for (uint32_t i = 0; i < no_bytes; ++i)
	{
		dest.push_back(data & 0xff);
		data >>= 8;
	}
}

// ************************************************************************************
void StoreUInt(uchar_t *dest, uint64_t data, uint32_t no_bytes)
{
	switch(no_bytes)
	{
	case 8:
		*dest++ = (data >> 56);
	case 7:
		*dest++ = (data >> 48) & 0xFF;
	case 6:
		*dest++ = (data >> 40) & 0xFF;
	case 5:
		*dest++ = (data >> 32) & 0xFF;
	case 4:
		*dest++ = (data >> 24) & 0xFF;
	case 3:
		*dest++ = (data >> 16) & 0xFF;
	case 2:
		*dest++ = (data >>  8) & 0xFF;
	case 1:
		*dest++ = (data      ) & 0xFF;
	}
}

// ************************************************************************************
void StoreUIntLSB(uchar_t *dest, uint64_t data, uint32_t no_bytes)
{
	for (uint32_t i = 0; i < no_bytes; ++i)
	{
		*dest++ = data & 0xff;
		data >>= 8;
	}
}

// ************************************************************************************
void StoreInt32LSB(uchar_t *dest, int32_t data, uint32_t no_bytes)
{
	for (uint32_t i = 0; i < no_bytes; ++i)
	{
		*dest++ = data & 0xff;
		data >>= 8;
	}
}

// ************************************************************************************
void StoreFloat(uchar_t* dest, float data)
{
	uchar_t* mem = (uchar_t*) alloca(5);
	memcpy(mem, (void *) &data, 4);

	for (uint32_t i = 0; i < 4; ++i)
		*dest++ = mem[i];
}

// ************************************************************************************
void LoadUInt(uchar_t *dest, uint64_t &data, uint32_t no_bytes)
{
	data = 0;
	switch(no_bytes)
	{
	case 8:
		data = (data << 8) + *dest++;
	case 7:
		data = (data << 8) + *dest++;
	case 6:
		data = (data << 8) + *dest++;
	case 5:
		data = (data << 8) + *dest++;
	case 4:
		data = (data << 8) + *dest++;
	case 3:
		data = (data << 8) + *dest++;
	case 2:
		data = (data << 8) + *dest++;
	case 1:
		data = (data << 8) + *dest++;
	}
}

// ************************************************************************************
void LoadUInt2(uchar_t *dest, uint32_t &data)
{
	data = *dest++;
	data = (data << 8) + *dest++;
}

// ************************************************************************************
void LoadUInt4(uchar_t *dest, uint32_t &data)
{
	data = *dest++;
	data = (data << 8) + *dest++;
	data = (data << 8) + *dest++;
	data = (data << 8) + *dest++;
}

// ************************************************************************************
void IncrementUInt(uchar_t *dest, uint32_t no_bytes)
{
	for(int i = no_bytes-1; i >= 0; --i)
		if(dest[i] < 0xff)
		{
			++dest[i];
			break;
		}
		else
			dest[i] = 0;
}

// ************************************************************************************
// Integer logarithm of base 
uint32_t IntLog2(uint64_t x)
{
	uint32_t r;

	for(r = 0; x > 1; ++r)
		x /= 2;

	return r;
}

// ************************************************************************************
// Integer logarithm of base 4
uint32_t IntLog4(uint64_t x)
{
	uint32_t r;

	for(r = 0; x > 1; ++r)
		x /= 4;

	return r;
}

// ************************************************************************************
// Convert integer to string
string Int2String(uint64_t x)
{
	stringstream tmp;
	tmp << x;

	string str;
	tmp >> str;

	return str;
}

// ************************************************************************************
int Int2PChar(uint64_t x, uchar_t *str)
{
	return CNumericConversions::Int2PChar(x, str);
}

// ************************************************************************************
int SInt2PChar(int64_t x, uchar_t *str)
{
	int r = 0;

	if(x < 0)
	{
		*str++ = '-';
		r++;
		x = -x;
	}

	return CNumericConversions::Int2PChar(x, str) + r;
}

// ************************************************************************************
int SDouble2PChar(double x, uint32_t prec, uchar_t *str)
{
	int r = 0;

	if (x < 0)
	{
		*str++ = '-';
		r++;
		x = -x;
	}

	return CNumericConversions::Double2PChar(x, prec, str) + r;
}

// ************************************************************************************
string Int2StringFilled(uint64_t x, uint32_t n)
{
	string s = Int2String(x);

	while(s.size() < n)
		s = "0" + s;

	return s;
}

// ************************************************************************************
string FormatInt(uint64_t x)
{
	string r;
	char s[20];

	while(x)
	{
		long long tmp = x % 1000;
		x /= 1000;
		sprintf(s, "%lld", tmp);
		string ss(s);
		if(x)
			while(ss.size() < 3)
				ss = "0" + ss;
		if(r.size())
			r = ss + "," + r;
		else
			r = ss;
	}

	return r;
}

// ************************************************************************************
// Return string describing stage
string StageDesc(uint32_t stage_major, uint32_t stage_minor, bool sensitive_mode)
{
	char s[256];
	if(sensitive_mode)
		sprintf(s, "stage %ds:%d", stage_major, stage_minor);
	else
		sprintf(s, "stage %d:%d", stage_major, stage_minor);

	return string(s);
}

// ************************************************************************************
// Return stage id encoded as: top bits (major) and 5 lowest bits (minor)
string StageDesc(uint32_t stage_id)
{
	return StageDesc((stage_id >> 5) & 0x1f, stage_id & 0x1f, (bool) (stage_id >> 10));
}

// ************************************************************************************
// Normalize directory name
string NormalizeDirectory(string dir)
{
	if(dir.empty())
	{
		dir = "./";
		return dir;
	}

	for(uint32_t i = 0; i < dir.size(); ++i)
		if(dir[i] == '\\')
			dir[i] = '/';

	if(dir.back() != '/')
		dir.push_back('/');

	return dir;
}

// ************************************************************************************
uchar_t DNA2Bin(uchar_t x)
{
	switch (x)
	{
	case 'A':
		return sym_code_A;
	case 'C':
		return sym_code_C;
	case 'G':
		return sym_code_G;
	case 'T':
		return sym_code_T;
	}

	return sym_code_N_ref;
}

// ************************************************************************************
// Convert prefix to unsigned int
uint32_t Prefix2Int(string prefix)
{
	uint32_t res = 0;

	for(auto c: prefix)
		switch(c)
		{
		case 'A':
			res = (res << 2) + sym_code_A;
			break;
		case 'C':
			res = (res << 2) + sym_code_C;
			break;
		case 'G':
			res = (res << 2) + sym_code_G;
			break;
		case 'T':
			res = (res << 2) + sym_code_T;
			break;
		}

	return res;
}

// ************************************************************************************
// Check size of the file
int64_t FileSize(string file_name)
{
	FILE *f = fopen(file_name.c_str(), "rb");
	
	if (!f)
		return -1;

	my_fseek(f, 0, SEEK_END);
	int64_t size = my_ftell(f);

	fclose(f);

	return size;
}

// ************************************************************************************
// ************************************************************************************

// ************************************************************************************
CProgress::CProgress()
{
	max_value = 0;
}

// ************************************************************************************
CProgress::~CProgress()
{

}

// ************************************************************************************
void CProgress::Init(int64_t _max_value, bool _show_comment)
{
	unique_lock<mutex> lck(mtx);

	max_value = _max_value;
	current_value = 0;

	show_comment = _show_comment;
}

// ************************************************************************************
void CProgress::SetComment(string _comment)
{
	unique_lock<mutex> lck(mtx);

	comment = _comment;
}

// ************************************************************************************
void CProgress::Step(int64_t increment)
{
	unique_lock<mutex> lck(mtx);

	current_value += increment;

	show_progress();
}

// ************************************************************************************
void CProgress::show_progress(void)
{
	stringstream stext;

	stext.width(5);
	stext.precision(1);
	stext.setf(ios_base::fixed);

	stext << 100.0 * current_value / max_value;
	stext << "%";

	if (show_comment)
		stext << " " << comment;

	stext << "\r";

	string text = stext.str();

	if (text != prev_text)
	{
		cerr << text;
		fflush(stderr);
		prev_text = text;
	}
}

// EOF
