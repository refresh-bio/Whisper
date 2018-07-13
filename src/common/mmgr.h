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


#ifndef _MMGR_H
#define _MMGR_H

#include "../common/defs.h"
#include <memory>
#include <condition_variable>
#include <mutex>
#include <iostream>
#include <algorithm>

//#define DISABLE_MEMORY_MONITOR

using namespace std;

// ************************************************************************************
// CPtrPool
// ************************************************************************************
class CPtrPool
{
	typedef struct {
		uchar_t *ptr;
		uint64_t size;
		bool reserved;
	} pool_item_t;
	
	pool_item_t *pool;
	int pool_size;
	int no_reserved;

	int64_t alloc_mem;

	mutable mutex mtx;							// The mutex to synchronise on
	condition_variable cv;						// The condition to wait for

public:
	CPtrPool(int _pool_size)
	{
		alloc_mem = 0;

		pool_size = _pool_size;
		pool = new pool_item_t[pool_size];
		no_reserved = 0;

		for (int i = 0; i < pool_size; ++i)
		{
			pool[i].size = 1;
			pool[i].ptr = new uchar_t[pool[i].size];
			pool[i].reserved = false;

			alloc_mem += 1;
		}
	}

	~CPtrPool()
	{
		for (int i = 0; i < pool_size; ++i)
		{
			delete[] pool[i].ptr;

			alloc_mem -= pool[i].size;
//			cerr << "Ptr_pool destructor: " + to_string(alloc_mem >> 20) + "MB\n";
		}

		delete[] pool;
	}

	uchar_t* Allocate(uint64_t requested_size)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this] {return no_reserved < pool_size; });

		// Try to find matching block
		for (int i = 0; i < pool_size; ++i)
			if (!pool[i].reserved && pool[i].size >= requested_size)
			{
				pool[i].reserved = true;
				++no_reserved;
				return pool[i].ptr;
			}

		for (int i = pool_size - 1; i >= 0; --i)
			if (!pool[i].reserved)
			{
				delete[] pool[i].ptr;

				alloc_mem -= pool[i].size;
//				cerr << "Ptr_pool alloc- : " + to_string(alloc_mem >> 20) + "MB\n";

				// 10% overhead to avoid too many reallocations
				pool[i].size = (uint64_t)(requested_size * 1.1);

				pool[i].ptr = new uchar_t[pool[i].size];
				pool[i].reserved = true;
				++no_reserved;
				alloc_mem += pool[i].size;
//				cerr << "Ptr_pool alloc+ : " + to_string(alloc_mem >> 20) + "MB\n";

				auto p = pool[i].ptr;

				sort(pool, pool + pool_size, [](pool_item_t &x, pool_item_t &y) {return x.size < y.size; });

				return p;
			}

//		cerr << "Cannot allocate in ptr pool: " + to_string(no_reserved) + " : " + to_string(pool_size) + "\n";
//		exit(1);

		return nullptr;
	}

	void Release(uchar_t *p)
	{
		lock_guard<mutex> lck(mtx);
		
		for (int i = 0; i < pool_size; ++i)
			if (pool[i].ptr == p)
			{
				pool[i].reserved = false;

				if (no_reserved-- == pool_size)
					cv.notify_one();
				return;
			}
	
//		cerr << "Cannot find pointer\n";
//		exit(1);
	}
};


// ************************************************************************************
// CMemoryPool
// ************************************************************************************
template<typename T> class CMemoryPool 
{
	int64_t total_size;
	int64_t part_size;
	int64_t n_parts_total;
	int64_t n_parts_free;

	uchar_t *buffer, *raw_buffer;
	uint32_t *stack;

	mutable mutex mtx;							// The mutex to synchronise on
	condition_variable cv;						// The condition to wait for

public:

	int64_t getAllocatedParts() { return n_parts_total - n_parts_free; }

	// ************************************************************************************
	CMemoryPool(int64_t _total_size, int64_t _part_size)
	{
		raw_buffer = nullptr;
		buffer = nullptr;
		stack  = nullptr;
		Prepare(_total_size, _part_size);
	}

	// ************************************************************************************
	~CMemoryPool()
	{
		Release();
	}

	// ************************************************************************************
	void Prepare(int64_t _total_size, int64_t _part_size)
	{
		Release();

		n_parts_total = _total_size / _part_size;
		part_size     = (_part_size + 15) / 16 * 16;			// to allow mapping pointer to int*
		n_parts_free  = n_parts_total;
	
		total_size = n_parts_total * part_size;

		raw_buffer = new uchar_t[total_size+64];
		buffer     = raw_buffer;
		while(((uint64_t) buffer) % 64)
			buffer++;

		stack = new uint32_t[n_parts_total];
		for(uint32_t i = 0; i < n_parts_total; ++i)
			stack[i] = i;
	}

	// ************************************************************************************
	void Release(void)
	{
		if(raw_buffer)
			delete[] raw_buffer;
		raw_buffer = nullptr;
		buffer     = nullptr;

		if(stack)
			delete[] stack;
		stack = nullptr;
	}

	// ************************************************************************************
	int64_t GetPartSize()
	{
		return part_size;
	}

	// ************************************************************************************
	int64_t GetAvailableParts()
	{
//		lock_guard<mutex> lck(mtx);

		return n_parts_free;
	}

	// ************************************************************************************
	void Reserve(T* &part)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 0;});
			
		part = reinterpret_cast<T*>(buffer + stack[--n_parts_free]*part_size);
	}

	// ************************************************************************************
	void Reserve(T* &part1, T* &part2)
	{
		unique_lock<mutex> lck(mtx);
		cv.wait(lck, [this]{return n_parts_free > 1;});
			
		part1 = reinterpret_cast<T*>(buffer + stack[--n_parts_free]*part_size);
		part2 = reinterpret_cast<T*>(buffer + stack[--n_parts_free]*part_size);
	}

	// ************************************************************************************
	void Free(T* part)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32_t) ((((uchar_t*) part) - buffer) / part_size);
		if(n_parts_free <= 2)
			cv.notify_all();		
	}

	// ************************************************************************************
	void FreeConditional(T* part)
	{
		if (part == nullptr) {
			return;
		}

		Free(part);
	}

	// ************************************************************************************
	void Free(T* part1, T* part2)
	{
		lock_guard<mutex> lck(mtx);

		stack[n_parts_free++] = (uint32_t) ((((uchar_t*) part1) - buffer) / part_size);
		stack[n_parts_free++] = (uint32_t) ((((uchar_t*) part2) - buffer) / part_size);
		if(n_parts_free <= 2)
			cv.notify_all();		
	}
};


// ************************************************************************************
// CMemoryMonitor
// ************************************************************************************
class CMemoryMonitor 
{
	uint64_t max_memory;
	uint64_t memory_in_use;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_memory_full;				// The condition to wait for

	void show_status(int64_t n)
	{
		cerr << "Memory monitor: " << max_memory << "  " << memory_in_use << "  : " << n << "\n";
	}

public:
	// ************************************************************************************
	CMemoryMonitor(uint64_t _max_memory)
	{
#ifdef DISABLE_MEMORY_MONITOR
		return;
#endif

		lock_guard<mutex> lck(mtx);
		max_memory    = _max_memory;
		memory_in_use = 0;
	}

	// ************************************************************************************
	~CMemoryMonitor()
	{}

	// ************************************************************************************
	// Try to increase the amount of allocated memory
	void Increase(uint64_t n) 
	{
#ifdef DISABLE_MEMORY_MONITOR
		return;
#endif

		unique_lock<mutex> lck(mtx);

		cv_memory_full.wait(lck, [this, n]{return memory_in_use + n <= max_memory;});
		memory_in_use += n;
	}

	// ************************************************************************************
	// Try to increase the amount of allocated memory
	// If no memory is accupied so far, completes even if the amount of requested memory is higher than the limit!
	void ForceIncrease(uint64_t n) 
	{
#ifdef DISABLE_MEMORY_MONITOR
		return;
#endif

		unique_lock<mutex> lck(mtx);

		cv_memory_full.wait(lck, [this, n]{return memory_in_use + n <= max_memory || memory_in_use == 0;});
		memory_in_use += n;
	}

	// ************************************************************************************
	// Decrease the amount of allocated memory
	void Decrease(uint64_t n) 
	{
#ifdef DISABLE_MEMORY_MONITOR
		return;
#endif

		lock_guard<mutex> lck(mtx);

		memory_in_use -= n;
		cv_memory_full.notify_all();
	}

	// ************************************************************************************
	// Return an info about the limit and the amount of allocated memory
	void Info(uint64_t &_max_memory, uint64_t &_memory_in_use)
	{
#ifdef DISABLE_MEMORY_MONITOR
		_max_memory = 0;
		_memory_in_use = 0;
		return;
#endif

		lock_guard<mutex> lck(mtx);
	
		_max_memory    = max_memory;
		_memory_in_use = memory_in_use;
	}
};

#endif

// EOF
