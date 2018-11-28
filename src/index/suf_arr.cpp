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


#include "suf_arr.h"
#include "suf_arr_core.h"
#include "../common/utils.h"
#include "../common/timer.h"

#include <iostream>
#include <memory>
#include <algorithm>

// ************************************************************************************
// Construct suffix array for direct genome
bool suffix_array_construct_dir(const string index_name, const string dest_dir, const string temp_dir)
{
	cerr << "\n***** Stage 5. Suffix array construction (direct)\n";

	shared_ptr<CMapperFile> in(new CMapperFile(EXT_REF_SEQ_DIR, MARKER_REF_SEQ_DIR));
	if(!in->OpenRead(temp_dir + index_name))
	{
		cerr << "Cannot open reference sequence (direct)\n";
		return false;
	}
	int64_t size = in->GetSize();

	shared_ptr<CMapperFile> out(new CMapperFile(EXT_SA_DIR, MARKER_SA_DIR, sizeof(uint32_t)));
	if(!out->OpenWrite(dest_dir + index_name))
	{
		cerr << "Cannot create suffix array (direct)\n";
		return false;
	}

	suffix_array_construct(in, out, (uint32_t) size, false);

	return true;
}

// ************************************************************************************
// Construct suffix array for rev. comp. genome
bool suffix_array_construct_rc(const string index_name, const string dest_dir, const string temp_dir)
{
	cerr << "\n***** Stage 6. Suffix array construction (rev. comp.)\n";

	shared_ptr<CMapperFile> in(new CMapperFile(EXT_REF_SEQ_RC, MARKER_REF_SEQ_RC));
	if(!in->OpenRead(temp_dir + index_name))
	{
		cerr << "Cannot open reference sequence (rev. comp.)\n";
		return false;
	}
	int64_t size = in->GetSize();

	shared_ptr<CMapperFile> out(new CMapperFile(EXT_SA_RC, MARKER_SA_RC, sizeof(uint32_t)));
	if(!out->OpenWrite(dest_dir + index_name))
	{
		cerr << "Cannot create suffix array (rev. comp.)\n";
		return false;
	}

	suffix_array_construct(in, out, (uint32_t) size, true);
	
	return true;
}

// ************************************************************************************
// Construct suffix array
bool suffix_array_construct(shared_ptr<CMapperFile> in, shared_ptr<CMapperFile> out, uint32_t size, bool translate_Ns)
{
	cerr << "Reading...\n";
	++size;
	shared_ptr<uchar_t> ref_seq(new uchar_t[size]);
	in->Read(ref_seq.get(), size-1);
	ref_seq.get()[size-1] = 0;

	if(translate_Ns)
	{
		cerr << "Translation Ns...\n";
		for(int64_t i = 0; i < size-1; ++i)
			if(ref_seq.get()[i] == sym_code_N_ref)
				ref_seq.get()[i] = 0;
			else
				ref_seq.get()[i] += 1;
	}

	cerr << "Construction... ";
	shared_ptr<uint32_t> sa(new uint32_t[size]);
	SA_IS(ref_seq.get(), sa.get(), size, 255, sizeof(uchar_t), 0);
	cerr << "\n";

	if(translate_Ns)
	{
		cerr << "Reverse translation Ns...\n";
		for(int64_t i = 0; i < size-1; ++i)
			if(ref_seq.get()[i] == 0)
				ref_seq.get()[i] = sym_code_N_ref;
			else
				ref_seq.get()[i] -= 1;
	}

	cerr << "Cleanup";
	suffix_array_cleanup(ref_seq, sa, size, lut_long_prefix_len);

	cerr << "Writing...";
	uint32_t *sa_ptr = sa.get();
	while(size)
	{
		uint32_t to_write = 1 << 27;

		if(to_write > size)
			to_write = size;
		cerr << ".";	fflush(stderr);
		out->Write(sa_ptr, to_write);
		size -= to_write;
		sa_ptr += to_write;
	}
	cerr << "\n";

	return true;
}

// ************************************************************************************
// Remove from suffix array the indexes pointing to suffixes at which any of prefix_len initial symbols is N
bool suffix_array_cleanup(shared_ptr<uchar_t> text, shared_ptr<uint32_t> sa, uint32_t &size, uint32_t prefix_len)
{
	uint32_t idx = 0;
	uint32_t* sa_ptr = sa.get();
	uchar_t* text_ptr = text.get();

	CStopWatch watch;

	watch.StartTimer();

	// Copy 0th position which points sentinel
	sa_ptr[idx++] = sa_ptr[0];
	for(uint32_t i = 1; i < size; ++i)
	{
		bool is_N = false;
		for(uint32_t j = sa_ptr[i]; j < min(sa_ptr[i]+prefix_len, size); ++j)
			if(text_ptr[j] >= sym_code_N_ref)
				is_N = true;

		if(!is_N)
			sa_ptr[idx++] = sa_ptr[i];

		if((i & ((1 << 26) - 1)) == 0)
		{
			cerr << ".";
			fflush(stderr);
		}
	}
	watch.StopTimer();

	cerr << "\nCompacted from " << size << " to " << idx << " entries in " << watch.GetElapsedTime() << "s\n";

	size = idx;

	return true;
}

// ************************************************************************************
// Compute LUT for direct genome
bool lut_compute_dir(const string index_name, const string dest_dir, const string temp_dir)
{
	cerr << "\n***** Stage 7. LUT computation (direct)\n";

	shared_ptr<CMapperFile> in_text(new CMapperFile(EXT_REF_SEQ_DIR, MARKER_REF_SEQ_DIR));
	if(!in_text->OpenRead(temp_dir + index_name))
	{
		cerr << "Cannot open reference sequence (direct)\n";
		return false;
	}
	int64_t size_text = in_text->GetSize();

	shared_ptr<CMapperFile> in_sa(new CMapperFile(EXT_SA_DIR, MARKER_SA_DIR, sizeof(uint32_t)));
	if(!in_sa->OpenRead(dest_dir + index_name))
	{
		cerr << "Cannot open suffix array (direct)\n";
		return false;
	}
	int64_t size_sa = in_sa->GetSize();

	shared_ptr<CMapperFile> out_short(new CMapperFile(EXT_LUT_SHORT_DIR, MARKER_LUT_SHORT_DIR, sizeof(uint32_t)));
	if(!out_short->OpenWrite(dest_dir + index_name))
	{
		cerr << "Cannot create LUT short (direct)\n";
		return false;
	}

	shared_ptr<CMapperFile> out_long(new CMapperFile(EXT_LUT_LONG_DIR, MARKER_LUT_LONG_DIR, sizeof(uint32_t)));
	if(!out_long->OpenWrite(dest_dir + index_name))
	{
		cerr << "Cannot create LUT long (direct)\n";
		return false;
	}

	lut_compute(in_text, in_sa, out_short, out_long, (uint32_t) size_text, (uint32_t) size_sa);

	return true;
}

// ************************************************************************************
// Compute LUT for rev. comp. genome
bool lut_compute_rc(const string index_name, const string dest_dir, const string temp_dir)
{
	cerr << "\n***** Stage 8. LUT computation (rev. comp.)\n";

	shared_ptr<CMapperFile> in_text(new CMapperFile(EXT_REF_SEQ_RC, MARKER_REF_SEQ_RC));
	if(!in_text->OpenRead(temp_dir + index_name))
	{
		cerr << "Cannot open reference sequence (rev. comp.)\n";
		return false;
	}
	int64_t size_text = in_text->GetSize();

	shared_ptr<CMapperFile> in_sa(new CMapperFile(EXT_SA_RC, MARKER_SA_RC, sizeof(uint32_t)));
	if(!in_sa->OpenRead(dest_dir + index_name))
	{
		cerr << "Cannot open suffix array (rev. comp.)\n";
		return false;
	}
	int64_t size_sa = in_sa->GetSize();

	shared_ptr<CMapperFile> out_short(new CMapperFile(EXT_LUT_SHORT_RC, MARKER_LUT_SHORT_RC, sizeof(uint32_t)));
	if(!out_short->OpenWrite(dest_dir + index_name))
	{
		cerr << "Cannot create LUT short (rev. comp.)\n";
		return false;
	}

	shared_ptr<CMapperFile> out_long(new CMapperFile(EXT_LUT_LONG_RC, MARKER_LUT_LONG_RC, sizeof(uint32_t)));
	if(!out_long->OpenWrite(dest_dir + index_name))
	{
		cerr << "Cannot create LUT long (rev. comp.)\n";
		return false;
	}

	lut_compute(in_text, in_sa, out_short, out_long, (uint32_t) size_text, (uint32_t) size_sa);

	return true;
}

// ************************************************************************************
// Computation of LUT
bool lut_compute(shared_ptr<CMapperFile> in_text, shared_ptr<CMapperFile> in_sa, shared_ptr<CMapperFile> out_short, shared_ptr<CMapperFile> out_long, uint32_t size_text, uint32_t size_sa)
{
	shared_ptr<uchar_t> ref_seq(new uchar_t[size_text]);
	shared_ptr<uint32_t> sa(new uint32_t[size_sa]);

	cerr << "Reading reference sequence...\n";
	in_text->Read(ref_seq.get(), size_text);
	cerr << "Reading suffix array";
	uint32_t to_read;
	uint32_t size = size_sa;
	uint32_t *ptr = sa.get();
	
	while(size)
	{
		to_read = 1 << 27;
		if(to_read > size)
			to_read = size;
		in_sa->Read(ptr, to_read);
		size -= to_read;
		ptr += to_read;
		cerr << ".";
		fflush(stderr);
	}
	cerr << "\n";

	uint32_t size_lut_short = (1 << (2 * lut_short_prefix_len)) + 1;
	uint32_t size_lut_long  = (1 << (2 * lut_long_prefix_len)) + 1;
	shared_ptr<uint32_t> lut_short(new uint32_t[size_lut_short]);
	shared_ptr<uint32_t> lut_long(new uint32_t[size_lut_long]);

	cerr << "Computing LUT";
	uint32_t lut_value_prev = ~(0u);

	uchar_t *text = ref_seq.get();
	uint32_t *sa_ptr = sa.get();

	// Skip 0th position (sentinel)
	for(uint32_t i = 1; i < size_sa; ++i)
	{
		uint32_t lut_value = 0;
		for(uint32_t j = 0; j < lut_long_prefix_len; ++j)
			lut_value = (lut_value << 2) + text[sa_ptr[i]+j];

		while(lut_value != lut_value_prev)
		{
			lut_long.get()[++lut_value_prev] = i;
		}

		if((i & ((1 << 26) - 1)) == 0)
		{
			cerr << ".";
			fflush(stderr);
		}
	}
	while(lut_value_prev != size_lut_long-1)
		lut_long.get()[++lut_value_prev] = size_sa;

	uint32_t lut_scale = (size_lut_long - 1) / (size_lut_short - 1);
	for(uint32_t i = 0; i < size_lut_short; ++i)
		lut_short.get()[i] = lut_long.get()[i * lut_scale];
	cerr << "\n";

	cerr << "Writing LUT...\n";
	out_short->Write(lut_short.get(), size_lut_short);
	out_long->Write(lut_long.get(), size_lut_long);
	
	return true;
}

// EOF
