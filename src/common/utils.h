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


#ifndef _UTILS_H
#define _UTILS_H

#include "defs.h"
#include <functional>
#include <cstring>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <mutex>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include "../libs/zlib.h"
#include "../libs/libdeflate.h"

#undef min
#undef max

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.707106781186547524401
#endif

using namespace std;

//#define NORM(x, min_val, max_val)	((x) < (min_val) ? (min_val) : (x) > (max_val) ? (max_val) : (x))
#define BITS2BYTES(x)			(((x) + 7) / 8)
#define PACKED_READ_LEN(x)		(((x) + 1) / 2)
#define GET_SYMBOL(seq, pos)	(((seq)[(pos)/2] >> ((((pos)+1) % 2) * 4)) & 0xf)
#define ASSERT(cnd, msg)		// Assert(cnd, msg)

//#define MEMCMP(x, y, len)		A_memcmp(x, y, len)
#define MEMCMP(x, y, len)		memcmp(x, y, len)

// ************************************************************************************
template<typename T>
inline T NormalizeValue(T x, T min_value, T max_value)
{
	if (x < min_value)
		return min_value;
	else if (x > max_value)
		return max_value;
	return x;
}

// ************************************************************************************
inline void Assert(bool condition, const std::string& msg) 
{
	if (!condition) {
		std::cerr << msg << std::endl;
		std::flush(std::cerr);
		exit(1);
	}
}

// ************************************************************************************
inline double log_base(double x, double base) 
{
	return log10(x) / log10(base);
}

// Our pointer aligner - necessary as std::align works properly for (alignment > 8) only since C++17
void *ptr_align(void *ptr, size_t alignment);
void *alloc_aligned(void *&raw_ptr, size_t size, size_t alignment);

typedef uint64_t read_id_t;

// ************************************************************************************
void StoreUIntLSB(vector<uchar_t> &dest, uint64_t data, uint32_t no_bytes);
void StoreUInt(uchar_t *dest, uint64_t data, uint32_t no_bytes);
void StoreUIntLSB(uchar_t *dest, uint64_t data, uint32_t no_bytes);
void StoreInt32LSB(uchar_t *dest, int32_t data, uint32_t no_bytes);
void LoadUInt(uchar_t *dest, uint64_t &data, uint32_t no_bytes);
void LoadUInt2(uchar_t *dest, uint32_t &data);
void LoadUInt4(uchar_t *dest, uint32_t &data);
void IncrementUInt(uchar_t *dest, uint32_t no_bytes);
uint32_t IntLog2(uint64_t x);
uint32_t IntLog4(uint64_t x);
string Int2String(uint64_t x);
string Int2StringFilled(uint64_t x, uint32_t n);
string FormatInt(uint64_t x);
int Int2PChar(uint64_t x, uchar_t *str);
int SInt2PChar(int64_t x, uchar_t *str);
int SDouble2PChar(double x, uint32_t prec, uchar_t *str);

string NormalizeDirectory(string dir);
uint32_t Prefix2Int(string prefix);

string StageDesc(uint32_t stage_major, uint32_t stage_minor, bool sensitive_mode = false);
string StageDesc(uint32_t stage_id);

int64_t FileSize(string file_name);


// ************************************************************************************
// Wrapper for FILE structure. 
// Added: 
//   * verification of markers (at the beginning and end of the file)
//   * can be used as parameter of shared_ptr template
// ************************************************************************************
class CMapperFile
{
	FILE *file;

	enum class open_mode_t {not_set, closed, opened_for_read, opened_for_write, direct_stream};

	open_mode_t state;

	string name;
	string extension;
	string marker;
	uint32_t element_size;
	int64_t size;
	uint64_t marker_length;

	bool file_error(string error_code);
	inline int64_t seek(int64_t offset, int origin);
	inline int64_t tell();

public:
	CMapperFile();
	CMapperFile(string _extension, string _marker, uint32_t _element_size = 1);
	~CMapperFile();

	bool SetParams(string _extension, string _marker, uint32_t _element_size = 1);
	bool OpenRead(string _name);
	bool OpenWrite(string _name);
	bool OpenStream(FILE *stream);
	bool Close();
	bool Remove();

	int64_t GetSize();
	int64_t GetRawSize();
	int64_t Read(void* data, int64_t count);
	int64_t ReadString(string &data);
	int64_t WriteByte(uint32_t data);
	int64_t WriteUint32(uint32_t data);
	int64_t Write(void* data, int64_t count);
	int64_t Write(const void* data, int64_t count);
	int64_t SeekBegin();
	int64_t Seek(int64_t pos);
	bool Eof();
	void SetBufSize(uint64_t size);
};

// ************************************************************************************
class CNumericConversions 
{
	public:
	static uchar_t digits[100000*5];
	static uint64_t powers10[15];
	struct _si {
		_si()
		{
			for(int i = 0; i < 100000; ++i)
			{
				int dig = i;
		
				digits[i*5+4] = '0' + (dig % 10);
				dig /= 10;
				digits[i*5+3] = '0' + (dig % 10);
				dig /= 10;
				digits[i*5+2] = '0' + (dig % 10);
				dig /= 10;
				digits[i*5+1] = '0' + (dig % 10);
				dig /= 10;
				digits[i*5+0] = '0' + dig;
			}	

			powers10[0] = 1;
			for (int i = 1; i < 15; ++i)
				powers10[i] = 10 * powers10[i - 1];
		}
	} static _init;
	
	static int NDigits(uint64_t v)
	{
		return (v < 10000)
			? (v < 100 ? (v < 10 ? 1 : 2) : (v < 1000 ? 3 : 4))
			: (v < 1000000 ? (v < 100000 ? 5 : 6) : (v < 10000000 ? 7 : 8));
	}

	static int Int2PChar(uint64_t val, uchar_t *str)
	{
		if(val >= 1000000000000000ull)
		{
			uint64_t dig1 = val / 1000000000000000ull;
			val -= dig1 * 1000000000000000ull;
			uint64_t dig2 = val / 10000000000ull;
			val -= dig2 * 10000000000ull;
			uint64_t dig3 = val / 100000ull;
			uint64_t dig4 = val - dig3 * 100000ull;
		
			int ndig = NDigits(dig1);
			
			memcpy(str, digits+dig1*5+(5-ndig), ndig);
			memcpy(str+ndig, digits+dig2*5, 5);
			memcpy(str+ndig+5, digits+dig3*5, 5);
			memcpy(str+ndig+10, digits+dig4*5, 5);
		
			return ndig+15;
		}
		else if(val >= 10000000000ull)
		{
			uint64_t dig1 = val / 10000000000ull;
			val -= dig1 * 10000000000ull;
			uint64_t dig2 = val / 100000ull;
			uint64_t dig3 = val - dig2 * 100000ull;
		
			int ndig = NDigits(dig1);

			memcpy(str, digits+dig1*5+(5-ndig), ndig);
			memcpy(str+ndig, digits+dig2*5, 5);
			memcpy(str+ndig+5, digits+dig3*5, 5);

			return ndig+10;
		}
		else if(val >= 100000ull)
		{
			uint64_t dig1 = val / 100000ull;
			uint64_t dig2 = val - dig1 * 100000ull;
		
			int ndig = NDigits(dig1);

			memcpy(str, digits+dig1*5+(5-ndig), ndig);
			memcpy(str+ndig, digits+dig2*5, 5);

			return ndig+5;
		}
		else
		{
			int ndig = NDigits(val);

			memcpy(str, digits+val*5+(5-ndig), ndig);

			return ndig;
		}
	}

	static int Double2PChar(double val, uint32_t prec, uchar_t *str)
	{
		int64_t a = (int64_t) val;
		int64_t b = (int64_t) ((1.0 + (val - (double) a)) * powers10[prec] + 0.5);
		
		int r1 = Int2PChar(a, str);
		int r2 = Int2PChar(b, str + r1);
		str[r1] = '.';
		
		return r1 + r2;
	}
};

// ************************************************************************************
// Call I/O operations in serial mode
class CSerialProcessing
{
	bool call_in_serial_mode;
	mutex mtx;

public:
	CSerialProcessing(bool _call_in_serial_mode = true) : call_in_serial_mode(_call_in_serial_mode)
	{}

	void Call(function<void(void)> &&call)
	{
		if (call_in_serial_mode)
		{
			unique_lock<mutex> lck(mtx);
			call();
		}
		else
			call();
	}
};

// ************************************************************************************
class CProgress
{
	mutex mtx;
	int show_comment;

	int64_t max_value;
	int64_t current_value;
	string comment;

	string prev_text;

	void show_progress(void);

public:
	CProgress();
	~CProgress();

	void Init(int64_t _max_value, bool _show_comment);
	void SetComment(string _comment);

	void Step(int64_t increment);
};

// ************************************************************************************
class CGzipMember
{
#ifndef _DEBUG
	z_stream stream;
#endif
	int compression_level;

public:
	CGzipMember(int _compression_level = 1) : compression_level(_compression_level)
	{
#ifndef _DEBUG
		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
#endif
	}

	~CGzipMember()
	{
	}

	size_t Compress(uchar_t *src, size_t src_size, uchar_t *dest, size_t dest_size)
	{
#ifndef _DEBUG
		stream.avail_in = (uInt)src_size;
		stream.next_in = (Bytef *)src;
		stream.avail_out = (uInt)dest_size;
		stream.next_out = (Bytef *)dest;

		deflateInit2(&stream, compression_level, Z_DEFLATED, MAX_WBITS + 16, 8, Z_DEFAULT_STRATEGY);

		if (deflate(&stream, Z_FINISH) != Z_STREAM_END)
			return 0;

		deflateEnd(&stream);

		return stream.total_out;
#else
		return 0;
#endif
	}
};

// ************************************************************************************
class CBGZF
{
	int compression_level;

	libdeflate_compressor *compressor;

public:
	CBGZF(int _compression_level) : compression_level(_compression_level ? _compression_level : 1)
	{
		compressor = libdeflate_alloc_compressor(compression_level);
	}

	~CBGZF()
	{
		libdeflate_free_compressor(compressor);
	}

	size_t Compress(uchar_t *src, size_t src_size, uchar_t *dest, size_t dest_size)
	{
		auto orig_dest = dest;

		*dest++ = 31;			// gzip ID
		*dest++ = 139;
		*dest++ = 8;			// compression method
		*dest++ = 4;			// flags
		
		*dest++ = 0;			// time (uint32)
		*dest++ = 0;
		*dest++ = 0;
		*dest++ = 0;

		*dest++ = 0;			// extra flags
		*dest++ = 0xff;			// OS

		*dest++ = 6;			// extra length (uint16)
		*dest++ = 0;

		*dest++ = 'B';			// BGZF id
		*dest++ = 'C';	
		*dest++ = 2;			// subfield size (uint16)
		*dest++ = 0;

		auto bsize_ptr = dest;
		dest += 2;

		auto deflate_size = libdeflate_deflate_compress(compressor, src, src_size, dest, dest_size - (dest - orig_dest) - 8);
		dest += deflate_size;

		StoreUIntLSB(bsize_ptr, deflate_size + 6 + 19, 2);

		uint32_t crc32 = 0;
		crc32 = libdeflate_crc32(crc32, src, src_size);

		StoreUIntLSB(dest, crc32, 4);
		dest += 4;

		StoreUIntLSB(dest, src_size, 4);
		dest += 4;

		return dest - orig_dest;
	}
};

#endif

// EOF
