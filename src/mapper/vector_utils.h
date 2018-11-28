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


#pragma once
#define MAX_VECTOR_SIZE 256
#include "../libs/vectorclass.h"

#include <string>
#include <sstream>
#include <iomanip>

enum class instruction_set_t { none, sse, sse2, sse3, sse3s, sse41, sse42, avx, avx2 };


inline std::string printVec2uq(const Vec2uq& v) {
	std::ostringstream oss;
	oss << std::hex << std::setw(16) << v.extract(0) << ", " << std::setw(16) << v.extract(1);;
	return oss.str();
}

inline std::string printVec4uq(const Vec4uq& v) {
	std::ostringstream oss;
	oss << std::hex << std::setw(16) << v.extract(0) << ", " << std::setw(16) << v.extract(1) << ", "
		<< std::setw(16) << v.extract(2) << ", " << std::setw(16) << v.extract(3);
	return oss.str();
}


