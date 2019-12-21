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


#include "sam_writer.h"
#include "../libs/asmlib.h"
#include <cstdlib>

//#define DISABLE_IO

// ************************************************************************************
// CSamWriter
// ************************************************************************************
CSamWriter::CSamWriter(CParams *params, CObjects *objects, vector<string> &header_SAM, vector<pair<string, uint32_t>> &header_BAM)
{
	// !!! Open SAM files
	total_size_mapped = 0;
	verbosity_level = params->verbosity_level;
	max_size = params->sam_buffer_memory;
	q_sam_blocks = objects->q_sam_blocks;
	gzipped_SAM_level = params->gzipped_SAM_level;
	store_BAM = params->store_BAM;
	use_stdout = params->use_stdout;

	buffer_mapped = new uchar_t[max_size];

	gzip = new CGzipMember(gzipped_SAM_level);

	mp_sam_parts = objects->mp_sam_parts;
	serial_processing = objects->serial_processing;

	if (!serial_processing)
	{
		cerr << "No serial processing object!\n";
		exit(1);
	}

	string file_name_mapped = params->project_name;

	if (store_BAM)
	{
		file_name_mapped += ".bam";
	}
	else
	{
		file_name_mapped += ".sam";

		if (gzipped_SAM_level)
			file_name_mapped += ".gz";
	}

	writes_no = 0;

	if (use_stdout)
		file_mapped.OpenStream(stdout);
	else
		file_mapped.OpenWrite(file_name_mapped);

	if (store_BAM)
		write_BAM_header(header_SAM, header_BAM);
	else
		write_SAM_header(header_SAM);
}

// ************************************************************************************
CSamWriter::~CSamWriter()
{
	delete gzip;
	delete[] buffer_mapped;
}

// ************************************************************************************
void CSamWriter::operator()()
{
	sam_block_t block;

	while (!q_sam_blocks->IsCompleted())
	{
		if (q_sam_blocks->Pop(block))
		{
			if (block.type == sam_results_t::mapped || block.type == sam_results_t::undefined)
			{
				if (total_size_mapped + block.size > max_size)
				{
#ifndef DISABLE_IO
					serial_processing->Call([&] {
						file_mapped.Write(buffer_mapped, total_size_mapped);
					});
#endif
					if (verbosity_level >= 3)
						cerr << "SamWriter: " << ++writes_no << " of size " << total_size_mapped << "\n";
					total_size_mapped = 0;
				}

				A_memcpy(buffer_mapped + total_size_mapped, block.data, block.size);
				total_size_mapped += block.size;
			}

			mp_sam_parts->Free(block.data);
		}
	}

	complete();
}

// ************************************************************************************
void CSamWriter::write_SAM_header(vector<string> &header_SAM)
{
	if (gzipped_SAM_level > 0)
		gz_begin_member(file_mapped, total_size_mapped, buffer_mapped);

	for (auto &x : header_SAM)
		if (gzipped_SAM_level == 0)
			write(file_mapped, total_size_mapped, buffer_mapped, x);
		else
			gz_write(file_mapped, total_size_mapped, buffer_mapped, x);

	if (gzipped_SAM_level > 0)
		gz_end_member(file_mapped, total_size_mapped, buffer_mapped);
}

// ************************************************************************************
void CSamWriter::write_BAM_header(vector<string> &header_SAM, vector<pair<string, uint32_t>> &header_BAM)
{
	// Preparing BAM header in binary uncompressed form
	vector<uchar_t> raw_header;

	raw_header.push_back('B');
	raw_header.push_back('A');
	raw_header.push_back('M');
	raw_header.push_back(1);

	uint32_t meta_len = 0;
	for (auto &x : header_SAM)
		meta_len += (uint32_t) x.size();

	StoreUIntLSB(raw_header, meta_len, 4);
	for (auto &x : header_SAM)
		raw_header.insert(raw_header.end(), x.begin(), x.end());

	StoreUIntLSB(raw_header, header_BAM.size(), 4);

	for (auto &x : header_BAM)
	{
		StoreUIntLSB(raw_header, x.first.size() + 1, 4);
		raw_header.insert(raw_header.end(), x.first.begin(), x.first.end());
		raw_header.push_back(0);
		StoreUIntLSB(raw_header, x.second, 4);
	}

	uchar_t *h = raw_header.data();
	uchar_t *dest = new uchar_t[64 << 10];
	CBGZF bgzf(gzipped_SAM_level);

	for (uint32_t i = 0; i < raw_header.size(); i += (62 << 10))
	{
		auto src = h + i;
		uint32_t src_size = min((uint32_t)raw_header.size() - i, 62u << 10);

		auto dest_size = bgzf.Compress(src, src_size, dest, 64u << 10);
		file_mapped.Write(dest, dest_size);
	}

	delete[] dest;
}

// ************************************************************************************
// Directly write some string into SAM file
void CSamWriter::write(CMapperFile &file, uint64_t &total_size, uchar_t* buffer, string s)
{
	if (total_size + s.size() > max_size)
	{
		file.Write(buffer, total_size);
		total_size = 0;
		if (verbosity_level >= 3)
			cerr << "SamWriter: " << ++writes_no << "\n";
		if (s.size() > max_size)
			file.Write(s.c_str(), s.size());
		else
		{
			copy_n(s.c_str(), s.size(), buffer);
			total_size = s.size();
		}
	}
	else
	{
		copy_n(s.c_str(), s.size(), buffer + total_size);
		total_size += s.size();
	}
}

// ************************************************************************************
// Start member gzipped block
void CSamWriter::gz_begin_member(CMapperFile &file, uint64_t &total_size, uchar_t *buffer)
{
	if (total_size)
	{
		file.Write(buffer, total_size);
		total_size = 0;
	}

	gz_raw_buffer.clear();
}

// ************************************************************************************
// Close member gzipped block
void CSamWriter::gz_end_member(CMapperFile &file, uint64_t &total_size, uchar_t *buffer)
{
	if (!gz_raw_buffer.empty())
	{
		size_t out_size = gz_raw_buffer.size();
		out_size += (out_size >> 12) + (out_size >> 14) + (out_size >> 15) + 13 + 20;
		uchar_t *out_buffer = new uchar_t[out_size];

		size_t comp_size = gzip->Compress((uchar_t*)gz_raw_buffer.c_str(), gz_raw_buffer.size(), out_buffer, out_size);
		file.Write(out_buffer, comp_size);

		delete[] out_buffer;

		gz_raw_buffer.clear();
		gz_raw_buffer.shrink_to_fit();
	}
}

// ************************************************************************************
// Directly write some string into gzipped SAM file
void CSamWriter::gz_write(CMapperFile &file, uint64_t &total_size, uchar_t* buffer, string s)
{
#ifndef _DEBUG
	gz_raw_buffer += s;
#endif
}

// ************************************************************************************
void CSamWriter::complete()
{
#ifndef DISABLE_IO
	serial_processing->Call([&] {
		file_mapped.Write(buffer_mapped, total_size_mapped);
	});
#endif
	if (verbosity_level >= 3)
		cerr << "SamWriter: " << ++writes_no << "\n";
	if (verbosity_level >= 3)
		cerr << "SamWriter: " << ++writes_no << "\n";

	if (verbosity_level >= 3)
		cerr << "SamWriter complete\n";

	if (store_BAM)
	{
		// Store EOF Block
		file_mapped.WriteByte(0x1f);
		file_mapped.WriteByte(0x8b);
		file_mapped.WriteByte(0x08);
		file_mapped.WriteByte(0x04);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0xff);
		file_mapped.WriteByte(0x06);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x42);
		file_mapped.WriteByte(0x43);
		file_mapped.WriteByte(0x02);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x1b);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x03);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
		file_mapped.WriteByte(0x00);
	}

	file_mapped.Close();
}

// EOF
