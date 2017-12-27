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


#include "../common/stats.h"
#include <algorithm>
#include <numeric>

// ************************************************************************************
CRunningStats::CRunningStats()
{
}

// ************************************************************************************
CRunningStats::~CRunningStats()
{
}

// ************************************************************************************
bool CRunningStats::Register(uint32_t id, string description, running_stats_t type)
{
#ifndef COLLECT_STATS
	return true;
#endif

	unique_lock<mutex> lck(mtx);

	if(descriptions.find(id) != descriptions.end())
		return false;

	descriptions[id]		      = make_pair(description, type);
	inv_descriptions[description] = make_pair(id, type);

	return true;
}

// ************************************************************************************
string CRunningStats::format_int(uint64_t x)
{
	string r;
	char s[20];

	while(x)
	{
		int64_t tmp = x % 1000;
		x /= 1000;
		sprintf(s, "%lld", (long long int) tmp);
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
bool CRunningStats::AddTotals(uint32_t id, int64_t value)
{
#ifndef COLLECT_STATS
	return true;
#endif

	unique_lock<mutex> lck(mtx);
	auto p = descriptions.find(id);

	if(p == descriptions.end())
		return false;

	stat_int_totals[p->second.first] += value;

	return true;
}

// ************************************************************************************
bool CRunningStats::AddTotals(uint32_t id, double value)
{
#ifndef COLLECT_STATS
	return true;
#endif

	unique_lock<mutex> lck(mtx);
	auto p = descriptions.find(id);

	if(p == descriptions.end())
		return false;

	stat_float_totals[p->second.first] += value;

	return true;
}

// ************************************************************************************
bool CRunningStats::AddValues(uint32_t id, int64_t value)
{
#ifndef COLLECT_STATS
	return true;
#endif

	unique_lock<mutex> lck(mtx);
	auto p = descriptions.find(id);

	if(p == descriptions.end())
		return false;

	stat_int_values[p->second.first].push_back(value);

	return true;
}

// ************************************************************************************
bool CRunningStats::AddValues(uint32_t id, double value)
{
#ifndef COLLECT_STATS
	return true;
#endif

	unique_lock<mutex> lck(mtx);
	auto p = descriptions.find(id);

	if(p == descriptions.end())
		return false;

	stat_float_values[p->second.first].push_back(value);

	return true;
}

// ************************************************************************************
void CRunningStats::CreateReport(ostream &stream)
{
#ifndef COLLECT_STATS
	return;
#endif

	unique_lock<mutex> lck(mtx);

	stream << "****************************************************\n";
	stream << "Totals:\n";
	stream << "---------------------------------------------------\n";

	for(auto &p : stat_int_totals)
		stream << p.first << ": " << format_int(p.second) << "\n";
	for(auto &p : stat_float_totals)
		stream << p.first << ": " << p.second << "\n";
	stream << "\n";
	
	stream << "****************************************************\n";
	stream << "Averages:\n";
	for(auto &p : stat_int_values)
	{
		if(inv_descriptions[p.first].second != running_stats_t::averages)
			continue;
		stream << "---------------------------------------------------\n";
		stream << p.first << ":\n";
		sort(p.second.begin(), p.second.end());
		
		int64_t total = accumulate(p.second.begin(), p.second.end(), (int64_t) 0);
		stream << "  total :  " << format_int(total) << "\n";
		stream << "  number:  " << format_int(p.second.size()) << "\n";
		stream << "  average: " << format_int(total / p.second.size()) << "\n";
		stream << "  median:  " << format_int(p.second[p.second.size()/2]) << "\n";
		stream << "  min:     " << format_int(p.second.front()) << "\n";
		stream << "  max:     " << format_int(p.second.back()) << "\n";
	}

	for(auto &p : stat_float_values)
	{
		if(inv_descriptions[p.first].second != running_stats_t::averages)
			continue;
		stream << "---------------------------------------------------\n";
		stream << p.first << ":\n";
		sort(p.second.begin(), p.second.end());
		
		double total = accumulate(p.second.begin(), p.second.end(), (double) 0.0);
		stream << "  total :  " << total << "\n";
		stream << "  number:  " << format_int(p.second.size()) << "\n";
		stream << "  average: " << (double) total / p.second.size() << "\n";
		stream << "  median:  " << p.second[p.second.size()/2] << "\n";
		stream << "  min:     " << p.second.front() << "\n";
		stream << "  max:     " << p.second.back() << "\n";
	}
	stream << "\n";

	stream << "****************************************************\n";
	stream << "Lists:\n";
	for(auto &p : stat_int_values)
	{
		if(inv_descriptions[p.first].second != running_stats_t::lists)
			continue;
		stream << "---------------------------------------------------\n";
		stream << p.first << ":\n";
		sort(p.second.begin(), p.second.end());
		
		int64_t total = accumulate(p.second.begin(), p.second.end(), (int64_t) 0);
		stream << "  total :  " << format_int(total) << "\n";
		stream << "  number:  " << format_int(p.second.size()) << "\n";
		stream << "  average: " << format_int(total / p.second.size()) << "\n";
		stream << "  median:  " << format_int(p.second[p.second.size()/2]) << "\n";
		stream << "  min:     " << format_int(p.second.front()) << "\n";
		stream << "  max:     " << format_int(p.second.back()) << "\n";
		stream << "  values:\n";
		for(auto &q : p.second)
			stream << "    " << format_int(q) << "\n";
	}

	for(auto &p : stat_float_values)
	{
		if(inv_descriptions[p.first].second != running_stats_t::lists)
			continue;
		stream << "---------------------------------------------------\n";
		stream << p.first << ":\n";
		sort(p.second.begin(), p.second.end());
		
		double total = accumulate(p.second.begin(), p.second.end(), (double) 0.0);
		stream << "  total :  " << total << "\n";
		stream << "  number:  " << format_int(p.second.size()) << "\n";
		stream << "  average: " << (double) total / p.second.size() << "\n";
		stream << "  median:  " << p.second[p.second.size()/2] << "\n";
		stream << "  min:     " << p.second.front() << "\n";
		stream << "  max:     " << p.second.back() << "\n";
		stream << "  values:\n";
		for(auto &q : p.second)
			stream << "    " << q << "\n";
	}

	stream << "\n";
}

// EOF
