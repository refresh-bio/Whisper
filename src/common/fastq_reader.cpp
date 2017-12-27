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


#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>

#include "../common/defs.h"
#include "../common/timer.h"
#include "../common/fastq_reader.h"
#include "../libs/asmlib.h"


// ************************************************************************************
// CGZFile - wrapper for controled reading of gzipped FASTQ files
// ************************************************************************************

// ************************************************************************************
bool CGZFile::Open(string file_name)
{
	if (in)
		Close();

	in = fopen(file_name.c_str(), "rb");

	if (!in)
		return false;

	if(serial_processing)
		serial_processing->Call([&]
		{
			my_fseek(in, 0, SEEK_END);
			in_size = my_ftell(in);
			my_fseek(in, 0, SEEK_SET);
		});
	else
	{
		my_fseek(in, 0, SEEK_END);
		in_size = my_ftell(in);
		my_fseek(in, 0, SEEK_SET);
	}

	// Init stream structure
#ifndef _DEBUG
	stream.zalloc = Z_NULL;
	stream.zfree = Z_NULL;
	stream.opaque = Z_NULL;
	stream.avail_in = 0;
	stream.next_in = Z_NULL;
	if (inflateInit2(&stream, 31) != Z_OK)
	{
		cerr << "Error while reading gz file\n";
		exit(1);
	}
	stream.avail_in = 0;
	stream.next_in = buffer;
#endif

	in_buffer_size = 0;
	is_eof = in_size == 0;

	return true;
}

// ************************************************************************************
bool CGZFile::Close()
{
	if (in)
	{
		fclose(in);
		in = nullptr;
	}

	return true;
}

// ************************************************************************************
bool CGZFile::Eof()
{
	if (!in)
		return true;
#ifndef _DEBUG
	return is_eof && (stream.avail_in == 0);
#endif

	return is_eof;
}

// ************************************************************************************
bool CGZFile::Read(uchar_t *ptr, uint32_t to_read, uint32_t &readed, uint32_t &raw_readed)
{
#ifndef _DEBUG
	if (!in)
		return false;

	if (Eof())
		return false;

	stream.next_out = ptr;
	stream.avail_out = to_read;

	raw_readed = 0;

	int ret;
	do
	{
		if (!stream.avail_in)
		{
			if(serial_processing)
				serial_processing->Call([&] {
					in_buffer_size = fread(buffer, 1, buffer_size, in);
				});
			else
				in_buffer_size = fread(buffer, 1, buffer_size, in);

			raw_readed += (uint32_t) in_buffer_size;

			stream.avail_in = (uint32_t) in_buffer_size;
			stream.next_in = buffer;
		}

		ret = inflate(&stream, Z_NO_FLUSH);

		switch (ret)
		{
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR;     /* and fall through */
		case Z_DATA_ERROR:
			cerr << "Some error while reading gzip file\n";
			exit(1);
		case Z_MEM_ERROR:
			inflateEnd(&stream);
			return false;
		}

		if (ret == Z_STREAM_END)
		{
			bool multistream = stream.avail_in || !feof(in);
			if (!multistream)
			{
				inflateEnd(&stream);
				is_eof = true;
				break;
			}
			else //multiple streams in one file
			{
				if (inflateReset(&stream) != Z_OK)
				{
					cerr << "Error while reading gzip file\n";
					exit(1);
				}
			}
		}
	} while (stream.avail_out);

	readed = to_read - stream.avail_out;
#endif

	return true;
}


// ************************************************************************************
// CFastqReader	- reader class
// ************************************************************************************

uint32_t CFastqReader::OVERHEAD_SIZE = 1 << 16;

// ************************************************************************************
// Constructor of FASTA/FASTQ reader
// Parameters:
//    * _mm - pointer to memory monitor (to check the memory limits)
CFastqReader::CFastqReader(bool _is_fasta, uint32_t _part_size, CObjects *_objects)
{
	is_fasta        = _is_fasta;
	part_size       = _part_size;

	internal_buffer = new uchar_t[part_size];			// !!! Uwzglêdniæ w bilansie pamiêci

	// Input file mode (default: uncompressed)
	mode      = mode_t::m_plain;

	// Pointers to input files in various formats (uncompressed, gzip-compressed, bzip2-compressed)
	in		  = nullptr;
#ifndef _DEBUG
	in_gzip   = nullptr;
	in_bzip2  = nullptr;
	bzerror   = BZ_OK;
	in_gz_file = nullptr;

	gzip_buffer_size  = 1 << 26;
	bzip2_buffer_size = 1 << 24;
#endif

	serial_processing = _objects->serial_processing;
	progress = _objects->progress;

	instruction_set_t instruction_set;
	int x = instrset_detect();
	if (x >= 0 && x <= 8)
		instruction_set = (instruction_set_t)x;
	else if (x < 0)
		instruction_set = instruction_set_t::none;
	else
		instruction_set = instruction_set_t::avx2;

	switch (instruction_set) {
	case instruction_set_t::none:
	case instruction_set_t::sse:
		throw new std::runtime_error("SSE2 extensions required!");
	case instruction_set_t::sse2:
	case instruction_set_t::sse3:
	case instruction_set_t::sse3s:
	case instruction_set_t::sse41:
	case instruction_set_t::sse42:
		ptr_CountEOLs = CountEOLs128<instruction_set_t::sse2>;
		break;
	case instruction_set_t::avx:
		ptr_CountEOLs = CountEOLs128<instruction_set_t::avx>;
		break;
	case instruction_set_t::avx2:
		ptr_CountEOLs = CountEOLs256;
		break;
	}
}

// ************************************************************************************
// Destructor - close the files
CFastqReader::~CFastqReader()
{
	close_files();

	if(internal_buffer)
		delete[] internal_buffer;

#ifndef _DEBUG
	if (in_gz_file)
		delete in_gz_file;
#endif
}

// ************************************************************************************
void CFastqReader::close_files()
{
	if(mode == mode_t::m_plain)
	{
		if(in)
			fclose(in);
		in = nullptr;
	}
#ifndef _DEBUG
	else if(mode == mode_t::m_gzip)
	{
		if (in_gz_file)
		{
			delete in_gz_file;
			in_gz_file = nullptr;
		}
	}
	else if(mode == mode_t::m_bzip2)
	{
		if(in)
		{
			BZ2_bzReadClose(&bzerror, in_bzip2);
			fclose(in);
		}
		in = nullptr;
	}
#endif
}

// ************************************************************************************
// 
void CFastqReader::Restart()
{
	close_files();

	part_filled = 0;
	mode        = mode_t::m_plain;
#ifndef _DEBUG
	bzerror     = BZ_OK;
#endif
}

// ************************************************************************************
// Open the file
bool CFastqReader::Open(string _input_file_name)
{
#ifndef _DEBUG
//	if(in || in_gzip || in_bzip2)
	if (in || in_gz_file || in_bzip2)
#else
	if(in)
#endif
		return false;

	input_file_name = _input_file_name;

	// Set mode according to the extension of the file name
	if(input_file_name.size() > 3 && string(input_file_name.end()-3, input_file_name.end()) == ".gz")
		mode = mode_t::m_gzip;
	else if(input_file_name.size() > 4 && string(input_file_name.end()-4, input_file_name.end()) == ".bz2")
		mode = mode_t::m_bzip2;
	else
		mode = mode_t::m_plain;

	// Uncompressed file
	if(mode == mode_t::m_plain)	
	{
		if((in = fopen(input_file_name.c_str(), "rb")) == nullptr)
			return false;
	}
#ifndef _DEBUG
	// Gzip-compressed file
	else if(mode == mode_t::m_gzip)
	{
		in_gz_file = new CGZFile(gzip_buffer_size);
		if (!in_gz_file->Open(input_file_name))
			return false;
	}
	// Bzip2-compressed file
	else if(mode == mode_t::m_bzip2)
	{
		in = fopen(input_file_name.c_str(), "rb");
		if(!in)
			return false;
		setvbuf(in, nullptr, _IOFBF, bzip2_buffer_size);
		if((in_bzip2 = BZ2_bzReadOpen(&bzerror, in, 0, 0, nullptr, 0)) == nullptr)
		{
			fclose(in);
			return false;
		}
	}
#endif

	part_filled = 0;

	return true;
}

// ************************************************************************************
void CFastqReader::count_reads(uchar_t *part, uint32_t filled, uint32_t &no_reads, uint32_t &last_read_start_pos)
{
	// Count LF
	auto n_eols = (*ptr_CountEOLs)(part, filled);

	// Calculate no. of. records
	int lines_per_read = (is_fasta ? 2 : 4);
	no_reads = (uint32_t) (n_eols / lines_per_read);

	// Find last read start pos
	int eols_to_find = n_eols % lines_per_read + 1;

	uint32_t i = filled-1;
	while (i)
	{
		if (part[i] == 0xA)
		{
			--eols_to_find;
			if (!eols_to_find)
				break;
		}
		--i;
	}

	if (i)
		last_read_start_pos = i + 1;
	else
		last_read_start_pos = 0;

	return;

	// Naive way
	uchar_t id_symbol;
	if (is_fasta)
		id_symbol = '>';
	else
		id_symbol = '@';

	bool was_EOL = true;
	no_reads = 0;
	last_read_start_pos = 0;

	int line_no = 0;

	for (uint32_t i = 0; i < filled; ++i)
	{
		if (was_EOL)
		{
			if (part[i] == id_symbol && line_no % 2 == 0)
			{
				++no_reads;
				was_EOL = false;
				last_read_start_pos = i;
			}
		}

		if (part[i] == '\r' || part[i] == '\n')
		{
			if (!was_EOL)
				++line_no;
			was_EOL = true;
		}
		else
			was_EOL = false;
	}
}

// ************************************************************************************
uint32_t CFastqReader::find_nth_read(uchar_t *part, uint32_t size, uint32_t no_reads)
{
	uchar_t id_symbol;
	if (is_fasta)
		id_symbol = '>';
	else
		id_symbol = '@';

	bool was_EOL = true;
	uint32_t local_no_reads = 0;
	last_read_start_pos = 0;

	int line_no = 0;

	for (uint32_t i = 0; i < size; ++i)
	{
		if (was_EOL)
		{
			if (part[i] == id_symbol && line_no % 2 == 0)
			{
				++local_no_reads;
				was_EOL = false;
				last_read_start_pos = i;

				if (local_no_reads == no_reads+1)
					return i;
			}
		}

		if (part[i] == '\r' || part[i] == '\n')
		{
			if (!was_EOL)
				++line_no;
			was_EOL = true;
		}
		else
			was_EOL = false;
	}

	return size;
}

// ************************************************************************************
// Read a part of the file
bool CFastqReader::GetPartInfo(uchar_t *&part, uint32_t &size, uint32_t &no_reads)
{
#ifndef _DEBUG
	if (!in && !in_gz_file && !in_bzip2)
#else
	if (!in)
#endif
	{
		no_reads = 0;
		size = 0;
		return false;
	}

	if (Eof())
	{
		no_reads = 0;
		size = 0;
		return false;
	}

	uint32_t readed;
	uint32_t raw_readed;

	// Copy from internal buffer
//	copy(internal_buffer, internal_buffer+part_filled, part);
	A_memcpy(part, internal_buffer, part_filled);

	// Read data
	if (mode == mode_t::m_plain)
		raw_readed = readed = (uint32_t)fread(part + part_filled, 1, part_size - part_filled, in);
#ifndef _DEBUG
	else if (mode == mode_t::m_gzip)
		in_gz_file->Read(part + part_filled, (int)part_size - part_filled, readed, raw_readed);
	else if(mode == mode_t::m_bzip2)
		readed = BZ2_bzRead(&bzerror, in_bzip2, part+part_filled, (int) part_size-part_filled);
#endif
	else
		readed = 0;				// Never should be here

	if(progress)
		progress->Step(raw_readed);

	total_filled = part_filled + readed;

	if(Eof())
	{
		size = total_filled;
		part_filled = 0;

		count_reads(part, total_filled, no_reads, last_read_start_pos);
		part_no_reads = no_reads;
		last_read_start_pos = total_filled;

		return true;
	}
	
	count_reads(part, total_filled, no_reads, last_read_start_pos);
	part_no_reads = no_reads;
	
	size = last_read_start_pos;
	
	return true;
}

// ************************************************************************************
bool CFastqReader::GetPart(uchar_t *part, uint32_t &size, uint32_t no_reads)
{
	if (no_reads < part_no_reads)
		last_read_start_pos = find_nth_read(part, size, no_reads);

	A_memcpy(internal_buffer, part + last_read_start_pos, total_filled - last_read_start_pos);

	part_filled = total_filled - last_read_start_pos;

	size = last_read_start_pos;

	return size != 0;
}

// ************************************************************************************
// Skip to next EOL from the current position in a buffer
bool CFastqReader::SkipNextEOL(uchar_t *part, uint32_t &pos, uint32_t max_pos)
{
	uint32_t i;

	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' && part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;

	pos = i+1;

	return true;
}

// ************************************************************************************
// Check whether there is an EOF
bool CFastqReader::Eof()
{
	if(mode == mode_t::m_plain)
		return feof(in) != 0;
#ifndef _DEBUG
	else if(mode == mode_t::m_gzip)
		return in_gz_file->Eof();
	else if(mode == mode_t::m_bzip2)
		return bzerror == BZ_STREAM_END;
#endif

	return true;
}


// ************************************************************************************
// CPrePostBaseFastqReader - FastqReader thread class for preprocessing stage (base class)
// ************************************************************************************

// ************************************************************************************
CPrePostBaseFastqReader::CPrePostBaseFastqReader(bool _is_fasta, CParams *params, CObjects *objects, CIDStore *_id_store)
{
	fqr1 = new CFastqReader(_is_fasta, (uint32_t) params->block_size, objects);
	fqr2 = new CFastqReader(_is_fasta, (uint32_t) params->block_size, objects);

	mem_pool     = objects->mp_fastq_blocks;
	q_file_names = objects->q_file_names;
	id_store	 = _id_store;

	running_stats     = objects->running_stats;

	verbosity_level   = params->verbosity_level;
}

// ************************************************************************************
CPrePostBaseFastqReader::~CPrePostBaseFastqReader()
{
	if(fqr1)
		delete fqr1;
	if(fqr2)
		delete fqr2;
}

// ************************************************************************************
void CPrePostBaseFastqReader::operator()()
{
	file_name_no_t file_name_no;

	CThreadWatch thr_watch;
	thr_watch.StartTimer();

	while(!q_file_names->IsCompleted())
	{
		if(q_file_names->Pop(file_name_no))
		{
			fqr1->Restart();
			fqr2->Restart();
			bool opened = true;
			bool single_file = file_name_no.file_name2.empty();

			if(!(opened = fqr1->Open(file_name_no.file_name1)))
			{
				cerr << "Error: Cannot open file " << file_name_no.file_name1 << "\n";
				continue;
			}
			if(single_file)
				load_single_file(file_name_no);
			else
			{
				if(!opened || !fqr2->Open(file_name_no.file_name2))
				{
					cerr << "Error: Cannot open file " << file_name_no.file_name2 << "\n";
					continue;
				}
				load_pair_files_thr(file_name_no);
			}
		}
	}

	end_of_thread();

	thr_watch.StopTimer();

	store_stats(thr_watch.GetElapsedTime());
}

// ************************************************************************************
void CPrePostBaseFastqReader::load_single_file(file_name_no_t file_name_no)
{
	uchar_t *part;
	uint32_t part_filled;
	read_id_t id = empty_read_id;
	uint32_t part_no = 0;

	mem_pool->Reserve(part);
	if(verbosity_level >= 3)
		cout << "FastqReader: " << file_name_no.file_name1 << "  (free blocks: " << mem_pool->GetAvailableParts() << ")\n";
	
	uint32_t no_reads;

	while(fqr1->GetPartInfo(part, part_filled, no_reads))
	{
		fqr1->GetPart(part, part_filled, no_reads);

		process_single_block(file_name_no, part_no, id, part, part_filled);

		mem_pool->Reserve(part);
		if(verbosity_level >= 3)
			cout << "FastqReader: " << file_name_no.file_name1 << "  (free blocks: " << mem_pool->GetAvailableParts() << ")\n";

		++part_no;
	}

	mem_pool->Free(part);
	mark_eof(file_name_no, id);
}

// ************************************************************************************
void CPrePostBaseFastqReader::load_pair_files_thr(file_name_no_t file_name_no)
{
	mutex mtx;
	condition_variable cv1, cv2;

	uchar_t *part1, *part2;
	uint32_t part_filled1, part_filled2;

	read_id_t id1 = empty_read_id;
	read_id_t id2 = empty_read_id;
	uint32_t part_no = 0;

	uint32_t part1_no = 0;
	uint32_t part2_no = 0;

	uint32_t no_reads1, no_reads2;

	if(verbosity_level >= 3)
		cout << "FastqReader: " << file_name_no.file_name1 << " : " << file_name_no.file_name2 << "  (free blocks: " << mem_pool->GetAvailableParts() << ")\n";

	while (true)
	{
		bool p1_status, p2_status;
		
		// Reading Fastq parts
		mem_pool->Reserve(part1, part2);

		thread thr1a([&] {p1_status = fqr1->GetPartInfo(part1, part_filled1, no_reads1); });
		p2_status = fqr2->GetPartInfo(part2, part_filled2, no_reads2);
		thr1a.join();

		++part1_no;
		++part2_no;

		int no_reads = MIN(no_reads1, no_reads2);

		if ((p1_status ^ p2_status) || (no_reads1 != no_reads2 && no_reads == 0))
		{
			cout << "Error: Different number of reads in paired files!\n";
			exit(1);
		}

		if (!p1_status)
			break;

		thread thr1b([&] {fqr1->GetPart(part1, part_filled1, no_reads); });
		fqr2->GetPart(part2, part_filled2, no_reads);
		thr1b.join();

		process_pair_blocks(file_name_no, part_no, id1, part1, part_filled1, id2, part2, part_filled2);
	
		if (verbosity_level >= 3)
			cout << "FastqReader: " << file_name_no.file_name1 << " : " << file_name_no.file_name2 << "  (free blocks: " << mem_pool->GetAvailableParts() << ")\n";

		++part_no;
	}

	mem_pool->Free(part1, part2);
	mark_eof(file_name_no, id1);
}

// ************************************************************************************
void CPrePostBaseFastqReader::load_pair_files(file_name_no_t file_name_no)
{
	uchar_t *part1, *part2;
	uint32_t part_filled1, part_filled2;

	read_id_t id1 = empty_read_id;
	read_id_t id2 = empty_read_id;
	uint32_t part_no = 0;
	uint32_t no_reads1, no_reads2;

	// Reading Fastq parts
	mem_pool->Reserve(part1, part2);

	if(verbosity_level >= 3)
		cout << "FastqReader: " << file_name_no.file_name1 << " : " << file_name_no.file_name2 << "  (free blocks: " << mem_pool->GetAvailableParts() << ")\n";

	while(fqr1->GetPartInfo(part1, part_filled1, no_reads1) && fqr2->GetPartInfo(part2, part_filled2, no_reads2))
	{
		fqr1->GetPart(part1, part_filled1, no_reads1);
		fqr2->GetPart(part2, part_filled2, no_reads2);

		process_pair_blocks(file_name_no, part_no, id1, part1, part_filled1, id2, part2, part_filled2);

		mem_pool->Reserve(part1, part2);

		if(verbosity_level >= 3)
			cout << "FastqReader: " << file_name_no.file_name1 << " : " << file_name_no.file_name2 << "  (free blocks: " << mem_pool->GetAvailableParts() << ")\n";

		++part_no;
	}
	
	mem_pool->Free(part1, part2);
	mark_eof(file_name_no, id1);
}

// ************************************************************************************
// CPreFastqReader - FastqReader thread class for preprocessing stage 
// ************************************************************************************

// ************************************************************************************
void CPreFastqReader::process_single_block(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id, uchar_t* part, uint32_t part_filled)
{
	id = id_store->RegisterBlock(file_name_no.file_no, file_name_no.file_name1, part_no, id);
	q_blocks->Push(fastq_block_t(id, part, part_filled));
}

// ************************************************************************************
void CPreFastqReader::process_pair_blocks(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id1, uchar_t* part1, uint32_t part_filled1, 
	read_id_t &id2, uchar_t* part2, uint32_t part_filled2)
{
	id1 = id_store->RegisterBlock(file_name_no.file_no+0, file_name_no.file_name1, part_no, id1);
	id2 = id_store->RegisterBlock(file_name_no.file_no+1, file_name_no.file_name2, part_no, id2);

	q_blocks->Push(fastq_block_t(id1, part1, part_filled1));
	q_blocks->Push(fastq_block_t(id2, part2, part_filled2));
}

// ************************************************************************************
void CPreFastqReader::store_stats(double time)
{
	running_stats->AddValues(STAT_TIME_THR_FASTQ_READER, time);
}

// ************************************************************************************
void CPreFastqReader::end_of_thread()
{
	q_blocks->MarkCompleted();
}

// ************************************************************************************
// CPostFastqReader - FastqReader thread class for postprocessing stage
// ************************************************************************************

// ************************************************************************************
void CPostFastqReader::process_single_block(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id, uchar_t* part, uint32_t part_filled)
{
	read_id_t prev_id = id;

	id = id_store->GetBlockID(file_name_no.file_no, file_name_no.file_name1, part_no);

	// Check whether there is new results group
	if(is_new_group(prev_id, id))
		joiner_mgr->MarkAllBinBlocks((uint32_t) (prev_id >> (id_bits_local + id_bits_subgroup)));
	if(is_new_group(prev_id, id) || prev_id == empty_read_id)
		q_res_ids->Push((uint32_t) (id >> (id_bits_local + id_bits_subgroup)));

	joiner_mgr->PutFastqBlock(fastq_block_t(id, part, part_filled), true);
	send_bytes += part_filled;
}

// ************************************************************************************
void CPostFastqReader::process_pair_blocks(file_name_no_t file_name_no, uint32_t part_no, read_id_t &id1, uchar_t* part1, uint32_t part_filled1, 
	read_id_t &id2, uchar_t* part2, uint32_t part_filled2)
{
	read_id_t prev_id1 = id1;

	id1 = id_store->GetBlockID(file_name_no.file_no+0, file_name_no.file_name1, part_no);
	id2 = id_store->GetBlockID(file_name_no.file_no+1, file_name_no.file_name2, part_no);

	// Check whether there is new results group
	if(is_new_group(prev_id1, id1))
		joiner_mgr->MarkAllBinBlocks((uint32_t) (prev_id1 >> (id_bits_local + id_bits_subgroup)));
	if(is_new_group(prev_id1, id1) || prev_id1 == empty_read_id)
		q_res_ids->Push((uint32_t) (id1 >> (id_bits_local + id_bits_subgroup)));
	
	joiner_mgr->PutFastqBlock(fastq_block_t(id1, part1, part_filled1), false);
	joiner_mgr->PutFastqBlock(fastq_block_t(id2, part2, part_filled2), false);
	send_bytes += part_filled1;
	send_bytes += part_filled2;
}

// ************************************************************************************
void CPostFastqReader::store_stats(double time)
{
	running_stats->AddValues(STAT_TIME_THR_FASTQ_READER_PP, time);
}

// ************************************************************************************
void CPostFastqReader::mark_eof(file_name_no_t file_name_no, read_id_t id)
{
	joiner_mgr->MarkAllBinBlocks((uint32_t) (id >> (id_bits_local + id_bits_subgroup)));
}

// ************************************************************************************
bool CPostFastqReader::is_new_group(read_id_t prev_id, read_id_t new_id)
{
	if(prev_id == empty_read_id)
		return false;
	if(((prev_id >> id_bits_local) & sub_block_mask) == sub_block_mask)	// new group
		return true;
	else
		return false;
}

// ************************************************************************************
void CPostFastqReader::end_of_thread()
{
	if(verbosity_level >= 2)
		cout << "PostFastqReader end_of_thread\n";
	q_res_ids->MarkCompleted();
	joiner_mgr->MarkFastqReaderCompleted();
}

// EOF
