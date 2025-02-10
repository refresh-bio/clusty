// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "log.h"

#include <sstream>
#include <iomanip>

using namespace std;

const int Log::LEVEL_DEBUG = 0;
const int Log::LEVEL_VERBOSE = 1;
const int Log::LEVEL_NORMAL = 2;

// *****************************************************************************************
//
std::string Log::formatLargeNumber(uint64_t num, int minWidth) {
	std::string out = "";

	do {
		uint64_t part = num % 1000LL;
		num = num / 1000LL;

		if (num > 0) {
			std::ostringstream oss;
			oss << "," << std::setw(3) << std::setfill('0') << part;
			out = oss.str() + out;
		}
		else {
			out = std::to_string(part) + out;
		}

	} while (num > 0);

	int initialSpaces = minWidth - (int)out.length();

	if (initialSpaces > 0) {
		out = string(initialSpaces, ' ') + out;
	}

	return out;
}
