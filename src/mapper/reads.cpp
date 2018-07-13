// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.1
// Date    : 2018-07-10
// License : GNU GPL 3
// *******************************************************************************************


#include "reads.h"
#include "../common/utils.h"
#include <algorithm>
#include <utility>
#include "../libs/asmlib.h"

using namespace std;

// ************************************************************************************
// Reads Collector
// ************************************************************************************
CReadsCollector::CReadsCollector(CParams *params, CObjects *objects)
{
	mp_bins_write = objects->mp_bins_write;
	q_bins_write  = objects->q_bins_write;
	
	no_bins    = params->no_bins;

	constant_read_len     = params->constant_read_len;
	max_reads_compression = params->max_reads_compression;

	part_size = (uint32_t) mp_bins_write->GetPartSize();
	read_id_bytes = BITS2BYTES(bin_total_bits);
	running_stats = objects->running_stats;

	prepare_buffers();
}

// ************************************************************************************
CReadsCollector::~CReadsCollector()
{
	if(buffers)
	{
		for(uint32_t i = 0; i < no_bins; ++i)
			if(buffers[i])
				mp_bins_write->Free(buffers[i]);
		delete[] buffers;
	}
	if(buffer_sizes)
		delete[] buffer_sizes;
	if(read_counters)
		delete[] read_counters;
}

// ************************************************************************************
void CReadsCollector::prepare_buffers()
{
	buffers       = new uchar_t*[no_bins];
	buffer_sizes  = new uint32_t[no_bins];
	read_counters = new uint32_t[no_bins];

	fill(buffer_sizes , buffer_sizes+no_bins , 0);
	fill(read_counters, read_counters+no_bins, 0);

	for(uint32_t i = 0; i < no_bins; ++i)
		mp_bins_write->Reserve(buffers[i]);
}

// ************************************************************************************
// Put read to bin_id bin
bool CReadsCollector::Push(uint32_t bin_id, read_id_t id, uchar_t *data, uint32_t size, uint32_t best_error, bool is_packed)
{
	if(buffer_sizes[bin_id] + size/2 + 32 >= part_size)
	{
		q_bins_write->Push(reads_bin_t(bin_id, bin_id, buffers[bin_id], (uint64_t) buffer_sizes[bin_id], (uint64_t) buffer_sizes[bin_id], (uint64_t) read_counters[bin_id]));
		mp_bins_write->Reserve(buffers[bin_id]);
		buffer_sizes[bin_id]  = 0;
		read_counters[bin_id] = 0;
	}

	// Insert read into buffer
	if(is_packed)
		buffer_sizes[bin_id] += copy_read(buffers[bin_id]+buffer_sizes[bin_id], id, data, size, best_error);
	else
		buffer_sizes[bin_id] += pack_read(buffers[bin_id]+buffer_sizes[bin_id], id, data, size, best_error);

	++read_counters[bin_id];

	return true;
}

// ************************************************************************************
// Store all remaining parts of bins
bool CReadsCollector::Complete()
{
	for(uint32_t i = 0; i < no_bins; ++i)
	{
		if(buffer_sizes[i])
			q_bins_write->Push(reads_bin_t(i, i, buffers[i], (uint64_t) buffer_sizes[i], (uint64_t) buffer_sizes[i], (uint64_t) read_counters[i]));
		else
			// If bin is empty, it is not pushed to writing queue, so it must be dealocated
			mp_bins_write->Free(buffers[i]);

		buffer_sizes[i] = 0;
		read_counters[i] = 0;
		buffers[i] = nullptr;
	}

	return true;
}

// ************************************************************************************
// Format:
//   * id         - (default: 5B)
//   * raw_len    - read length: (2B)
//   * best_error - (default: 1B)
//   * data       - ceil(raw_len/2)
// Read is packed
inline uint32_t CReadsCollector::pack_read(uchar_t* dest, read_id_t id, uchar_t *data, uint32_t size, uint32_t best_error)
{
	uint32_t offset = 0;

	StoreUInt(dest, id, read_id_bytes);					// read id
	offset += read_id_bytes;
	if(!constant_read_len)
	{
		StoreUInt(dest+offset, size, 2);					// read len (raw)
		offset += 2;
	}
	StoreUInt(dest+offset, best_error, 1);				// best error
	offset++;

	if(max_reads_compression)							// 3 symbols in 1 byte
	{
		uint32_t buffer = 0;
		uint32_t buffer_bits = 0;
		uint32_t c;

		for(uint32_t i = 0; i < size/3; ++i)
		{
			c = data[i*3]*25 + data[i*3+1]*5 + data[i*3+2];
			buffer = (buffer << 7) + c;
			buffer_bits += 7;
			if(buffer_bits >= 8)
			{
				dest[offset++] = (buffer >> (buffer_bits - 8)) & 0xff;
				buffer_bits -= 8;
			}
		}

		if(size % 3 == 1)
		{
			buffer = (buffer << 3) + data[size-1];
			buffer_bits += 3;
		}
		else if(size % 3 == 2)
		{
			buffer = (buffer << 5) + data[size-2]*5 + data[size-1];
			buffer_bits += 5;
		}

		if(buffer_bits >= 8)
		{
			dest[offset++] = (buffer >> (buffer_bits - 8)) & 0xff;
			buffer_bits -= 8;
		}
		if(buffer_bits)
			dest[offset++] = (buffer << (8 - buffer_bits)) & 0xff;
	}
	else												// 2 symbols in 1 byte
	{
		for(uint32_t i = 0; i < size/2; ++i)
			dest[offset++] = (data[i*2] << 4) + data[i*2+1];

		if(size & 1)
			dest[offset++] = data[size-1] << 4;
	}

	return offset;
}

// ************************************************************************************
// Format:
//   * id         - (default: 5B)
//   * raw_len    - read length: (2B)
//   * best_error - (default: 1B)
//   * data       - ceil(raw_len/2)
// Read is just copied
inline uint32_t CReadsCollector::copy_read(uchar_t* dest, read_id_t id, uchar_t *data, uint32_t size, uint32_t best_error)
{
	uint32_t offset = 0;

	StoreUInt(dest, id, read_id_bytes);					// read id
	offset += read_id_bytes;
	if(!constant_read_len)
	{
		StoreUInt(dest+offset, size, 2);				// read len (raw)
		offset += 2;
	}
	StoreUInt(dest+offset, best_error, 1);				// best error
	offset++;

	if(max_reads_compression)							// 3 symbols in 1 byte
	{
		uint32_t buffer = 0;
		uint32_t buffer_bits = 0;
		uint32_t c;

		for(uint32_t i = 0; i < size/3; ++i)
		{
			c = GET_SYMBOL(data, i*3)*25 + GET_SYMBOL(data, i*3+1)*5 + GET_SYMBOL(data, i*3+2);
			buffer = (buffer << 7) + c;
			buffer_bits += 7;
			if(buffer_bits >= 8)
			{
				dest[offset++] = (buffer >> (buffer_bits - 8)) & 0xff;
				buffer_bits -= 8;
			}
		}

		if(size % 3 == 1)
		{
			buffer = (buffer << 3) + GET_SYMBOL(data, size-1);
			buffer_bits += 3;
		}
		else if(size % 3 == 2)
		{
			buffer = (buffer << 5) + GET_SYMBOL(data, size-2)*5 + GET_SYMBOL(data, size-1);
			buffer_bits += 5;
		}

		if(buffer_bits >= 8)
		{
			dest[offset++] = (buffer >> (buffer_bits - 8)) & 0xff;
			buffer_bits -= 8;
		}
		if(buffer_bits)
			dest[offset++] = (buffer << (8 - buffer_bits)) & 0xff;
	}
	else
	{
		for(uint32_t i = 0; i < size/2; ++i)
			dest[offset++] = data[i];

		if(size & 1)
			dest[offset++] = data[size/2];
	}

	return offset;
}

// ************************************************************************************
// CReadsDeliverer
// ************************************************************************************

// ************************************************************************************
CReadsDeliverer::CReadsDeliverer(CParams *params, uint32_t _read_id_len, uint32_t _start_pos, uint32_t _end_pos, CRunningStats *_running_stats,
	CMemoryPool<uchar_t> *_mp_reads, uint32_t _stage_id)
{
	is_valid      = false;
	is_sorted     = false;
	is_compressed = false;
	bin.data      = nullptr;

	read_id_len = _read_id_len;
	start_pos    = _start_pos;
	end_pos      = _end_pos;

	min_read_len          = params->min_read_len;
	max_read_len		  = params->max_read_len;

	constant_read_len     = params->constant_read_len;
	max_reads_compression = params->max_reads_compression;

	running_stats = _running_stats;
	stage_id      = _stage_id;

	mp_reads = _mp_reads;

	mp_reads->Reserve(buf_data_sft);
	
	if(max_reads_compression)
		prepare_lut();

	sorted_ptrs = nullptr;
}

// ************************************************************************************
CReadsDeliverer::~CReadsDeliverer()
{
	mp_reads->Free(buf_data_sft);

	if(sorted_ptrs)
		delete[] sorted_ptrs;
}

// ************************************************************************************
bool CReadsDeliverer::SetBin(reads_bin_t &_bin)
{
	bin       = _bin;
	is_valid  = true;
	is_sorted = false;

	is_compressed = max_reads_compression;

	if(sorted_ptrs)
		delete[] sorted_ptrs;
	sorted_ptrs = nullptr;

	return true;
}

// ************************************************************************************
uint64_t CReadsDeliverer::GetNoReads()
{
	if(!is_valid)
		return 0;

	return bin.count;
}

// ************************************************************************************
bool CReadsDeliverer::Start()
{
	if(!is_valid)
		return false;

	if(is_compressed)
		decompress_reads();

	if(!is_sorted)
		sort_reads();

	pos = 0;

	return true;
}

// ************************************************************************************
// Get next read from bin
bool CReadsDeliverer::Pop(read_id_t &id, uchar_t* &data, uchar_t* &data_sft, uint32_t &size, uint32_t &best_error)
{
	if(!is_sorted)
		return false;

	if(pos >= bin.count)
		return false;

	// Read id
	uint64_t tmp, raw_length;
	uchar_t *ptr = sorted_ptrs[pos];
	LoadUInt(ptr, tmp, read_id_len);
	ptr += read_id_len;
	id = tmp;

	// Read length (raw)
	if(constant_read_len)
		raw_length = min_read_len;
	else
	{
		LoadUInt(ptr, raw_length, 2);
		ptr += 2;
	}
	uint32_t packed_length = (uint32_t) PACKED_READ_LEN(raw_length);

	// best error
	best_error = *ptr++;

	// Packed sequence
	data = ptr;
	size = (uint32_t) raw_length;
	ptr += packed_length;

	++pos;

	// Shift packed sequence by half of a byte
	data_sft = buf_data_sft;

	data_sft[0] = data[0] >> 4;
	if((size & 1) == 0)
	{
		for(uint32_t i = 1; i < size/2; ++i)
			data_sft[i] = ((data[i-1] & 0xf) << 4) + (data[i] >> 4);
		data_sft[size/2] = (data[size/2-1] & 0xf) << 4;
	}
	else
	{
		for(uint32_t i = 1; i <= size/2; ++i)
			data_sft[i] = ((data[i-1] & 0xf) << 4) + (data[i] >> 4);
	}
	
	return true;
}

// ************************************************************************************
void CReadsDeliverer::prepare_lut()
{
	for(uint32_t i = 0; i < 125; ++i)
		decompress_lut[i] = ((i / 25) << 8) + (((i / 5) % 5) << 4) + (i % 5);
}

// ************************************************************************************
void CReadsDeliverer::decompress_reads()
{
	uint32_t i, j;
	uint64_t raw_len;

	uchar_t *curr_in_ptr  = bin.data + bin.size - bin.raw_size;
	uchar_t *curr_out_ptr = bin.data;

	for(i = 0; i < bin.count; ++i)
	{
		// Copy read id
		for(j = 0; j < read_id_len; ++j)
			*curr_out_ptr++ = *curr_in_ptr++;

		if(constant_read_len)
			raw_len = min_read_len;
		else
		{
			LoadUInt(curr_in_ptr, raw_len, 2);			// Load read len and copy this len
			*curr_out_ptr++ = *curr_in_ptr++;
			*curr_out_ptr++ = *curr_in_ptr++;
		}

		*curr_out_ptr++ = *curr_in_ptr++;				// copy best error value

		// Decompress read
		uint32_t in_buffer       = 0;
		uint32_t in_buffer_bits  = 0;
		uint32_t out_buffer      = 0;
		uint32_t out_buffer_bits = 0;

		for(j = 0; j < raw_len / 3; ++j)
		{
			if(in_buffer_bits < 7)
			{
				in_buffer = (in_buffer << 8) + *curr_in_ptr++;
				in_buffer_bits += 8;
			}
			uint32_t c = (in_buffer >> (in_buffer_bits - 7)) & 0x7f;
			in_buffer_bits -= 7;

			out_buffer = (out_buffer << 12) + decompress_lut[c];
			out_buffer_bits += 12;

			while(out_buffer_bits >= 8)
			{
				*curr_out_ptr++ = (out_buffer >> (out_buffer_bits - 8)) & 0xff;
				out_buffer_bits -= 8;
			}
		}

		if(raw_len % 3 == 1)
		{
			if(in_buffer_bits < 3)
			{
				in_buffer = (in_buffer << 8) + *curr_in_ptr++;
				in_buffer_bits += 8;
			}
			in_buffer >>= in_buffer_bits - 3;
			out_buffer = (out_buffer << 4) + (in_buffer & 0x07);
			out_buffer_bits += 4;
		}
		else if(raw_len % 3 == 2)
		{
			if(in_buffer_bits < 5)
			{
				in_buffer = (in_buffer << 8) + *curr_in_ptr++;
				in_buffer_bits += 8;
			}
			in_buffer >>= in_buffer_bits - 5;
			uint32_t cc = in_buffer & 0x1f;
			out_buffer = (out_buffer << 8) + ((cc / 5) << 4) + (cc % 5);
			out_buffer_bits += 8;
		}

		while(out_buffer_bits >= 8)
		{
			*curr_out_ptr++ = (out_buffer >> (out_buffer_bits - 8)) & 0xff;
			out_buffer_bits -= 8;
		}
		if(out_buffer_bits == 4)
			*curr_out_ptr++ = (out_buffer & 0x0f) << 4;
	}

	is_compressed = false;
}

// ************************************************************************************
void CReadsDeliverer::sort_reads()
{
	uchar_t** ptrs = new uchar_t*[bin.count];			// pointers to reads
	uchar_t* curr_ptr = bin.data;
	uint64_t read_id;
	uint32_t raw_len;
	uint32_t packed_len;
	uint32_t data_offset = read_id_len + (constant_read_len ? 0 : 2) + 1;
	uint32_t sorting_part_len;
	uint32_t start_offset = data_offset + start_pos / 2;

	if(start_pos % 2 == 0)
		sorting_part_len = (end_pos - start_pos + 2) / 2; 
	else
		sorting_part_len = (end_pos - start_pos + 1) / 2; 

	// Find pointers to the reads
	for(uint32_t i = 0; i < bin.count; ++i)
	{
		ptrs[i] = curr_ptr;
		LoadUInt(curr_ptr, read_id, read_id_len);
		curr_ptr += read_id_len;
		if(constant_read_len)
			raw_len = min_read_len;
		else
		{
			LoadUInt2(curr_ptr, raw_len);
			curr_ptr += 2;
		}
		packed_len = (uint32_t) PACKED_READ_LEN(raw_len);
		curr_ptr++;						// best error value
		curr_ptr += packed_len;
	}

	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	// Sort reads of contant lengths 
	if(constant_read_len || (max_read_len == min_read_len))
	{
		if(start_pos == 0)				// substage x:0
			sort(ptrs, ptrs+bin.count, [&](uchar_t* p, uchar_t* q){
				return MEMCMP(p+start_offset, q+start_offset, packed_len) < 0;
			});
		else if(start_pos % 2 == 0)		// data to sort start at byte boundary
			sort(ptrs, ptrs+bin.count, [&](uchar_t* p, uchar_t* q){
				int r = MEMCMP(p+start_offset, q+start_offset, sorting_part_len);
				if(!r)
					r = MEMCMP(p+data_offset, q+data_offset, packed_len);
				return r < 0;
			});
		else						// data to sort start inside a byte
			sort(ptrs, ptrs+bin.count, [&](uchar_t* p, uchar_t* q){
				if((p[start_offset] & 0xf) != (q[start_offset] & 0xf))
					return (p[start_offset] & 0xf) < (q[start_offset] & 0xf);
				
				int r = MEMCMP(p+start_offset+1, q+start_offset+1, sorting_part_len);
				if(!r)
					r = MEMCMP(p+data_offset, q+data_offset, packed_len);
				return r < 0;
		});
	}
	// Sort reads of different lengths 
	else
	{
		uint32_t min_packed_len = PACKED_READ_LEN(min_read_len);

		if(start_pos == 0)				// substage x:0
			sort(ptrs, ptrs+bin.count, [&](uchar_t* p, uchar_t* q){
				int r = MEMCMP(p+start_offset, q+start_offset, min_packed_len);

				if(r)
					return r < 0;

				uint32_t p_len, q_len;
				LoadUInt2(p+read_id_len, p_len);
				LoadUInt2(q+read_id_len, q_len);
				packed_len = (uint32_t) PACKED_READ_LEN(MIN(p_len, q_len));

				if(min_packed_len != packed_len)
				{
					r = MEMCMP(p+start_offset+min_packed_len, q+start_offset+min_packed_len, packed_len-min_packed_len);
					if(r)	
						return r < 0;
				}

				return p_len < q_len;
			});
		else if(start_pos % 2 == 0)		// data to sort start at byte boundary
			sort(ptrs, ptrs+bin.count, [&](uchar_t* p, uchar_t* q){
				int r = MEMCMP(p+start_offset, q+start_offset, sorting_part_len);
				if(r)
					return r < 0;

				r = MEMCMP(p+data_offset, q+data_offset, min_packed_len);
				if(r)
					return r < 0;

				uint32_t p_len, q_len;
				LoadUInt2(p+read_id_len, p_len);
				LoadUInt2(q+read_id_len, q_len);
				packed_len = (uint32_t) PACKED_READ_LEN(MIN(p_len, q_len));

				if(min_packed_len != packed_len)
				{
					r = MEMCMP(p+data_offset+min_packed_len, q+data_offset+min_packed_len, packed_len-min_packed_len);
					if(r)
						return r < 0;
				}

				return p_len < q_len;
			});
		else						// data to sort start inside a byte
			sort(ptrs, ptrs+bin.count, [&](uchar_t* p, uchar_t* q){
				if((p[start_offset] & 0xf) != (q[start_offset] & 0xf))
					return (p[start_offset] & 0xf) < (q[start_offset] & 0xf);
				
				int r = MEMCMP(p+start_offset+1, q+start_offset+1, sorting_part_len);
				if(r)
					return r < 0;

				r = MEMCMP(p+data_offset, q+data_offset, min_packed_len);
				if(r)
					return r < 0;

				uint32_t p_len, q_len;
				LoadUInt2(p+read_id_len, p_len);
				LoadUInt2(q+read_id_len, q_len);
				packed_len = (uint32_t) PACKED_READ_LEN(MIN(p_len, q_len));

				if(min_packed_len != packed_len)
				{
					r = MEMCMP(p+data_offset+min_packed_len, q+data_offset+min_packed_len, packed_len-min_packed_len);
					if(r)
						return r < 0;
				}

				return p_len < q_len;
		});
	}

	thr_watch.StopTimer();

	running_stats->AddValues(STAT_SORTING_TOTAL, thr_watch.GetElapsedTime());
	running_stats->AddValues(STAT_SORTING_BASE+stage_id, thr_watch.GetElapsedTime());

	sorted_ptrs = ptrs;

	is_sorted = true;
}

// ************************************************************************************
// CReadsSplitter
// ************************************************************************************
CReadsSplitter::CReadsSplitter(CParams *params, CObjects *objects, bool _is_fasta, uint32_t first_major_stage)
{
	mp_bins_write = objects->mp_bins_write;
	mp_reads      = objects->mp_reads;
	q_bins_write  = objects->q_bins_write;
	prefix_map    = params->prefix_map;
	no_bins       = prefix_map.back()+2;		// the last element of prefix map contains the last bin_id
	read_id_len   = BITS2BYTES(bin_total_bits);
	is_fasta      = _is_fasta;
 
	mask_low_quality_bases = params->mask_low_quality_bases;

	running_stats = objects->running_stats;

	bin_prefix = IntLog4(prefix_map.size());

	reads_collector = new CReadsCollector(params, objects);

	for(int i = 0; i < 256; ++i)
		codes[i] = sym_code_N_read;
	codes['A'] = codes['a'] = sym_code_A;
	codes['C'] = codes['c'] = sym_code_C;
	codes['G'] = codes['g'] = sym_code_G;
	codes['T'] = codes['t'] = sym_code_T;

	hist_read_len.resize(max_fastq_rec_length, 0);

	n_reads = 0;
}

// ************************************************************************************
CReadsSplitter::~CReadsSplitter()
{
	delete reads_collector;
}

// ************************************************************************************
bool CReadsSplitter::Process(uchar_t *_block, uint32_t _block_size, read_id_t id_range)
{
	read_id_t read_id = id_range;

	block      = _block;
	block_size = _block_size;
	block_pos  = 0;

	uchar_t *seq;
	uint32_t seq_size;

	mp_reads->Reserve(seq);

	while(get_seq(seq, seq_size))
	{
		// Calculation of bin id.
		// In a case of any N i prefix, the bin_id is the maximal value
		uint32_t bin_id = 0;
		if(seq_size < bin_prefix)
			bin_id = no_bins-1;
		else
		{
			uint32_t i;
			for (i = 0; i < seq_size / 2; ++i)
			{
				uchar_t sym = seq[i];
				if(sym == sym_code_N_read)
					break;
				if(i < bin_prefix)
					bin_id = (bin_id << 2) + sym;
			}

			if(i < seq_size / 2)
				bin_id = no_bins-1;
			else
				bin_id = prefix_map[bin_id];		// maps long prefix of a read onto bin id
		}

		reads_collector->Push(bin_id, read_id, seq, seq_size, error_unknown, false);
		read_id += 2;
		++n_reads;
		++hist_read_len[seq_size];
	}

	mp_reads->Free(seq);

	return true;
}

// ************************************************************************************
void CReadsSplitter::Complete()
{
	reads_collector->Complete();
}

// ************************************************************************************
inline bool CReadsSplitter::get_seq(uchar_t *seq, uint32_t &seq_size)
{
	uchar_t c;
	uint32_t pos = 0;

	if(is_fasta)
	{
		// Title
		if(block_pos >= block_size)
			return false;
		c = block[block_pos++];
		if(c != '>')
			return false;
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(block_pos >= block_size)
			return false;

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Sequence
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			seq[pos++] = codes[c];
		}
		seq_size = pos;

		if(block_pos >= block_size)
			return true;

		if(block[block_pos++] >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return true;
	}
	else
	{
		// Title
		if(block_pos >= block_size)
			return false;
		c = block[block_pos++];
		if(c != '@')
			return false;
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(block_pos >= block_size)
			return false;

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Sequence
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			seq[pos++] = codes[c];
		}
		if(block_pos >= block_size)
			return false;

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Plus
		c = block[block_pos++];
		if(block_pos >= block_size)
			return false;
		if(c != '+')
			return false;
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
		}
		if(block_pos >= block_size)
			return false;

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Quality
		if (mask_low_quality_bases >= 33)
		{
			if (block_pos + pos >= block_size)
				return false;

			for (uint32_t i = 0; i < pos; ++i)
			{
				c = block[block_pos++];
				if (c < mask_low_quality_bases)
					seq[i] = sym_code_N_read;
			}
		}
		else
		{
			block_pos += pos;
			if (block_pos >= block_size)
				return false;
		}

		c = block[block_pos++];
		seq_size = pos;

		if(block_pos >= block_size)
			return true;

		if(block[block_pos++] >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return true;
	}

	return (c == '\n' || c == '\r');
}

// ************************************************************************************
// CWReadsSplitter
// ************************************************************************************

// ************************************************************************************
CWReadsSplitter::CWReadsSplitter(CParams *params, CObjects *objects, bool _is_fasta, uint32_t first_major_stage)
{
	rs = new CReadsSplitter(params, objects, _is_fasta, first_major_stage);

	mp_blocks	   = objects->mp_fastq_blocks;
	q_fastq_blocks = objects->q_blocks;
	q_bins_write   = objects->q_bins_write;

	verbosity_level = params->verbosity_level;
}

// ************************************************************************************
CWReadsSplitter::~CWReadsSplitter()
{
	if(rs)
		delete rs;
}

// ************************************************************************************
void CWReadsSplitter::operator()()
{
	// Splitting parts
	while(!q_fastq_blocks->IsCompleted())
	{
		fastq_block_t fastq_block;
		if(q_fastq_blocks->Pop(fastq_block))
		{
			if(verbosity_level >= 2)
			{
				cerr << "Reads splitter: " << hex << fastq_block.id_range << dec << "\n";
			}
			rs->Process(fastq_block.data, fastq_block.size, fastq_block.id_range);
			mp_blocks->Free(fastq_block.data);
		}
	}
	rs->Complete();
	q_bins_write->MarkCompleted();

	rs->GetHistReadLen(hist_read_len);
}


// ************************************************************************************
// CReadsReader
// ************************************************************************************

// ************************************************************************************
CReadsReader::CReadsReader(CParams *params, CObjects *objects, bool _is_fasta)
{
	read_id_len   = BITS2BYTES(bin_total_bits);
	is_fasta      = _is_fasta;
 
	running_stats = objects->running_stats;

	n_reads = 0;

	mp_fastq_records = objects->mp_fastq_records;

	mp_fastq_records->Reserve(id);
	mp_fastq_records->Reserve(sequence);
	mp_fastq_records->Reserve(plus);
	mp_fastq_records->Reserve(quality);
}

// ************************************************************************************
CReadsReader::~CReadsReader()
{
	mp_fastq_records->Free(id);
	mp_fastq_records->Free(sequence);
	mp_fastq_records->Free(plus);
	mp_fastq_records->Free(quality);
}

// ************************************************************************************
bool CReadsReader::SetBlock(uchar_t *_block, uint64_t _block_size)
{
	block       = _block;
	block_size  = _block_size;
	block_pos   = 0;

	return block != nullptr;
}

// ************************************************************************************
bool CReadsReader::Pop(uchar_t* &_id, uchar_t* &_sequence, uchar_t* &_plus, uchar_t* &_quality,
	uint32_t &id_len, uint32_t &sequence_len, uint32_t &plus_len, uint32_t &quality_len)
{
	bool r = get_seq(id, sequence, plus, quality);

	_id       = id;
	_sequence = sequence;
	_plus     = plus;
	_quality  = quality;

	id_len       = (uint32_t) strlen((char*) id);
	sequence_len = (uint32_t) strlen((char*) sequence);
	plus_len     = (uint32_t) strlen((char*) plus);
	quality_len  = (uint32_t) strlen((char*) quality);

	return r;
}

// ************************************************************************************
void CReadsReader::Restart()
{
	block_pos = 0;
}

// ************************************************************************************
inline bool CReadsReader::get_seq(uchar_t *id, uchar_t *sequence, uchar_t *plus, uchar_t *quality)
{
	uchar_t c;
	uint32_t seq_len = 0;

	if(is_fasta)
	{
		// Title
		if(block_pos >= block_size)
			return false;
		c = block[block_pos++];
		if(c != '>')
			return false;
		*id++ = c;
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			*id++ = c;
		}
		*id = '\0';
		if(block_pos >= block_size)
			return false;

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Sequence
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			*sequence++ = c;
		}
		*sequence = '\0';

		if(block_pos >= block_size)
			return true;

		if(block[block_pos++] >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return true;

		*plus = '\0';
		*quality = '\0';
	}
	else
	{
		// Title
		if(block_pos >= block_size)
			return false;
		c = block[block_pos++];
		if(c != '@')
			return false;
		*id++ = c;
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			*id++ = c;
		}
		if(block_pos >= block_size)
			return false;
		*id = '\0';

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Sequence
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			*sequence++ = c;
			++seq_len;
		}
		if(block_pos >= block_size)
			return false;

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;
		*sequence = '\0';

		// Plus
		c = block[block_pos++];
		if(block_pos >= block_size)
			return false;
		if(c != '+')
			return false;
		*plus++ = c;
		for(; block_pos < block_size;)
		{
			c = block[block_pos++];
			if(c < 32)					// newliners
				break;
			*plus++ = c;
		}
		if(block_pos >= block_size)
			return false;
		*plus = '\0';

		c = block[block_pos++];
		if(c >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return false;

		// Quality
		for(uint32_t i = 0; i < seq_len; ++i)
			*quality++ = block[block_pos++];
		if(block_pos >= block_size)
			return false;
		*quality = '\0';
		c = block[block_pos++];

		if(block_pos >= block_size)
			return true;

		if(block[block_pos++] >= 32)
			block_pos--;
		else if(block_pos >= block_size)
			return true;
	}

	return (c == '\n' || c == '\r');
}

// EOF
