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


#ifndef _STATS_H
#define _STATS_H

#include "../common/defs.h"
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <mutex>
#include <utility>
#include <ostream>

using namespace std;

enum class running_stats_t {totals, averages, lists};

// ************************************************************************************
class CRunningStats 
{
	mutable mutex mtx;							// The mutex to synchronise on

	map<uint32_t, pair<string, running_stats_t>> descriptions;
	map<string, pair<uint32_t, running_stats_t>> inv_descriptions;

	map<string, int64_t> stat_int_totals;
	map<string, vector<int64_t>> stat_int_values;
	map<string, double> stat_float_totals;
	map<string, vector<double>> stat_float_values;

	string format_int(uint64_t x);

public:
	CRunningStats();
	~CRunningStats();

	bool Register(uint32_t id, string description, running_stats_t type);
	bool AddTotals(uint32_t id, int64_t value);
	bool AddTotals(uint32_t id, double value);
	bool AddValues(uint32_t id, int64_t value);
	bool AddValues(uint32_t id, double value);

	void CreateReport(ostream &stream);
};

#endif

// EOF
