#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek
*/

#include <memory>
#include <string>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <cctype>
#include <limits>
#include <cmath>


// ************************************************************************************
class Conversions
{
public:
	static char digits[100000 * 5];
	static uint64_t powers10[15];
	static double neg_powers10[15];
	struct _si {
		_si()
		{
			for (int i = 0; i < 100000; ++i)
			{
				int dig = i;

				digits[i * 5 + 4] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 3] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 2] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 1] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 0] = '0' + dig;
			}

			powers10[0] = 1;
			neg_powers10[0] = 1.0;
			for (int i = 1; i < 15; ++i) {
				powers10[i] = 10 * powers10[i - 1];
				neg_powers10[i] = 0.1 * neg_powers10[i - 1];
			}
		}
	} static _init;

	static int NDigits(uint64_t v)
	{
		return (v < 10000)
			? (v < 100 ? (v < 10 ? 1 : 2) : (v < 1000 ? 3 : 4))
			: (v < 1000000 ? (v < 100000 ? 5 : 6) : (v < 10000000 ? 7 : 8));
	}

	static int Int2PChar(uint64_t val, char *str)
	{
		if (val >= 1000000000000000ull)
		{
			uint64_t dig1 = val / 1000000000000000ull;
			val -= dig1 * 1000000000000000ull;
			uint64_t dig2 = val / 10000000000ull;
			val -= dig2 * 10000000000ull;
			uint64_t dig3 = val / 100000ull;
			uint64_t dig4 = val - dig3 * 100000ull;

			int ndig = NDigits(dig1);

			std::memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			std::memcpy(str + ndig, digits + dig2 * 5, 5);
			std::memcpy(str + ndig + 5, digits + dig3 * 5, 5);
			std::memcpy(str + ndig + 10, digits + dig4 * 5, 5);

			return ndig + 15;
		}
		else if (val >= 10000000000ull)
		{
			uint64_t dig1 = val / 10000000000ull;
			val -= dig1 * 10000000000ull;
			uint64_t dig2 = val / 100000ull;
			uint64_t dig3 = val - dig2 * 100000ull;

			int ndig = NDigits(dig1);

			std::memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			std::memcpy(str + ndig, digits + dig2 * 5, 5);
			std::memcpy(str + ndig + 5, digits + dig3 * 5, 5);

			return ndig + 10;
		}
		else if (val >= 100000ull)
		{
			uint64_t dig1 = val / 100000ull;
			uint64_t dig2 = val - dig1 * 100000ull;

			int ndig = NDigits(dig1);

			memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			memcpy(str + ndig, digits + dig2 * 5, 5);

			return ndig + 5;
		}
		else
		{
			int ndig = NDigits(val);

			memcpy(str, digits + val * 5 + (5 - ndig), ndig);

			return ndig;
		}
	}

	static int Double2PChar(double val, uint32_t prec, char *str)
	{
		int64_t a = (int64_t)val;
		int64_t b = (int64_t)((1.0 + (val - (double)a)) * powers10[prec] + 0.5);

		int r1 = Int2PChar(a, str);
		int r2 = Int2PChar(b, str + r1);
		str[r1] = '.';

		return r1 + r2;
	}

	static int String2PChar(const std::string& src, char* dst) {
		memcpy(dst, src.c_str(), src.size());
		return src.size();
	}


	static long int strtol(const char* str, char** endptr) {
		long int val = 0;
		char* p = (char*)str;
		bool is_negative = false;


		if (*p == '-')
		{
			is_negative = true;
			++p;
		}

		while (*p >= '0' && *p <= '9')
		{
			val = val * 10 + (*p++ - '0');
		}

		if (endptr)
			*endptr = p;

		return is_negative ? -val : val;
	}


	static double strtod(const char* str, char** endptr) {
		char* p = (char*)str;
		double r = 0.0;
		bool neg = false;
		if (*p == '-') {
			neg = true;
			++p;
		}
		while (*p >= '0' && *p <= '9') {
			r = (r * 10.0) + (*p - '0');
			++p;
		}
		if (*p == '.') {
			double f = 0.0;
			int n = 0;
			++p;
			while (*p >= '0' && *p <= '9') {
				f = (f * 10.0) + (*p - '0');
				++p;
				++n;
			}
			r += f * neg_powers10[n];
		}

		// optional exponential part
		if (*p == 'e' || *p == 'E') {
			
			++p;
			int exp = 0;
			bool exp_neg = false;

			if (*p == '-') {
				exp_neg = true;
				++p;
			}

			while (*p >= '0' && *p <= '9') {
				exp = exp * 10 + (*p++ - '0');
			}
			
			r *=  exp_neg ? neg_powers10[exp] : powers10[exp];
		}

		if (neg) {
			r = -r;
		}

		if (endptr) {
			*endptr = p;
		}

		return r;
	}
};


template <class T>
class FixedPoint {
public:
	static const T INVALID = std::numeric_limits<T>::min();

	static T fromString(const char* s, int& decimal_places)
	{
		bool neg = (*s == '-');
		T v = 0;
		if (neg) {
			++s;
		}

		// must start from digit 
		if (*s < '0' || *s > '9') {
			return INVALID;
		}

		decimal_places = 0;
		int decimal_inc = 0;

		for (int i = 0; ; ++i) {
			char c = *s;
			if (c >= '0' && c <= '9') {
				v = 10 * v + (c - '0');
				decimal_places += decimal_inc;
			}
			else if (c == '.' && decimal_places == 0) {
				decimal_inc = 1;
			}
			else if (c == 0 || c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r') {
				break;
			}
			else {
				return INVALID;
			}
			++s;
		}
		return (neg ? -v : v);
	}

	static int toString(T v, int decimal_places, char* dst) {
		
		char* p = dst;
		if (v < 0) {
			*p++ = '-';
			v = -v;
		}

		if (decimal_places > 0) {
			T before = v / Conversions::powers10[decimal_places];
			T after = v % Conversions::powers10[decimal_places];
			p += Conversions::Int2PChar(before, p);
			*p++ = '.';

			// FIXME: this works for at most 5 decimal places
			memcpy(p, Conversions::digits + after * 5 + (5 - decimal_places), decimal_places);
			p += decimal_places;		
		}
		else {
			p += Conversions::Int2PChar(v, p);
		}

		return p - dst;
	}

	static void alterPrecision(T& value, int cur_decimals, int ref_decimals) {

		// if current precision is less then desired
		if (cur_decimals < ref_decimals) {
			value *= Conversions::powers10[ref_decimals - cur_decimals];
		}

		// if current precision is larger then desired
		else if (cur_decimals > ref_decimals) {
			value /= Conversions::powers10[cur_decimals - ref_decimals];
		}
	}
};


// integral specialization
template <typename Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type* = nullptr>
int num2str(Integer val, char *out) {
	return Conversions::Int2PChar((uint64_t)val, out);
}

// floating point specialization
template <typename Floating, typename std::enable_if<std::is_floating_point<Floating>::value, int>::type* = nullptr>
int num2str(Floating val, char *out) {
	if (val == 0) {
		*out = '0';
		return 1;
	}
	return Conversions::Double2PChar((double)val, 6, out);
}

// pair specialization
template <typename T, typename U>
int num2str(const std::pair<T,U> val, char *out) {
	char* ptr = out;
	ptr += num2str(val.first, ptr);
	*ptr++ = ':';
	ptr += num2str(val.second, ptr);
	
	return ptr - out;
}

// collection specialization
template <typename T>
int num2str(const T* collection, size_t size, char delim, char* out) {
	char* ptr = out;
	for (size_t i = 0; i < size; ++i) {
		ptr += num2str(*collection++, ptr);
		*ptr++ = delim;
	}

	return ptr - out;
}

template <typename T>
int num2str_sparse(const T* collection, size_t size, char delim, char* out) {
	char* ptr = out;
	for (size_t i = 0; i < size; ++i, collection++) {
		if (*collection != 0) {
			ptr += num2str(i + 1, ptr);
			*ptr++ = ':';
			ptr += num2str(*collection, ptr);
			*ptr++ = delim;
		}
	}

	return ptr - out;
}
