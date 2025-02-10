#pragma once
// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <vector>
#include <ostream>
#include <string>

#include "conversion.h"


class InputBuffer {

public:
	const size_t size;
	char* const data;

	char* block_begin;
	char* block_end;

public:

	InputBuffer(size_t size) : size(size), data(new char[size + 1]), block_begin(data), block_end(nullptr) {}
	~InputBuffer() { delete[] data; }
};



// general variant
template <class T>
void value2buffer(const T& val, char* & buf) {  }


// explicit specializations
template <>
inline void value2buffer<std::string>(const std::string& val, char*& buf) {
	memcpy(buf, val.c_str(), val.size());
	buf += val.size();
}

template <>
inline void value2buffer<std::string_view>(const std::string_view& val, char*& buf) {
	memcpy(buf, val.data(), val.size());
	buf += val.size();
}

template <>
inline void value2buffer<std::string*>(std::string* const& val, char*& buf) {
	memcpy(buf, val->c_str(), val->size());
	buf += val->size();
}

template<>
inline void value2buffer<int>(const int& val, char*& buf) {
	buf += num2str(val, buf);
}

// generaltemplate
template <int First, int HowMany, class... Ts>
struct tuple2buffer {

	void operator() (const std::tuple<Ts...>& row, char separator, char*& p) const {

		value2buffer(std::get<First>(row), p);
		*p++ = separator;
		tuple2buffer<First + 1, HowMany - 1, Ts...>()(row, separator, p);

	}
};

// explicit specialization
template <int First, class... Ts>
struct tuple2buffer<First, 1, Ts...> {

	void operator() (const std::tuple<Ts...>& row, char separator, char*& p) const {
		value2buffer(std::get<First>(row), p);
	}
};



// buffered table writer
template <int NumCols, class... Ts>
void saveTableBuffered(
	std::ostream& ofs,
	const std::array<std::string, NumCols>& names,
	const std::vector<std::tuple<Ts...>>& values,
	char separator) {

	const int BUF_SIZE = 64LL << 20; // 64 MB output buffer
	int max_line_len = 1024;
	char* buffer = new char[BUF_SIZE];
	char* pend = buffer + BUF_SIZE;
	char* p = buffer;

	for (const auto& name : names) {
		value2buffer(name, p);
		*p++ = separator;
	}
	--p;
	*p++ = '\n';

	for (size_t i = 0; i < values.size(); ++i) {
		tuple2buffer<0, NumCols, Ts...>()(values[i], separator, p);
		*p++ = '\n';

		// save chunk to file
		if (pend - p < max_line_len) {
			ofs.write(buffer, p - buffer);
			p = buffer;
		}
	}

	ofs.write(buffer, p - buffer);

	delete[] buffer;

}


// buffered table writer
template <class T, class U>
void saveTableBuffered(
	std::ostream& ofs, 
	const std::string& name1,
	const std::vector<T>& column1, 
	const std::string name2,
	const std::vector<U>& column2, char separator) {

	const int BUF_SIZE = 64LL << 20; // 64 MB output buffer
	int max_line_len = 1024;
	char* buffer = new char[BUF_SIZE];
	char* pend = buffer + BUF_SIZE;
	char* p = buffer;

	value2buffer(name1, p);
	*p++ = separator;
	value2buffer(name2, p);
	*p++ = '\n';

	for (size_t i = 0; i < column1.size(); ++i) {
		value2buffer(column1[i], p);
		*p++ = separator;
		value2buffer(column2[i], p);
		*p++ = '\n';

		// save chunk to file
		if (pend - p < max_line_len) {
			ofs.write(buffer, p - buffer);
			p = buffer;
		}
	}

	ofs.write(buffer, p - buffer);

	delete[] buffer;

}

