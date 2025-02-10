#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <cstdint>
#include <string_view>


// *******************************************************************************************
class StringHasher {
	std::hash<char> hasher;
public:
	StringHasher() {}

	size_t operator()(const char* s) const {
		size_t hs = hasher(*s);
		++s;
		while (*s) {
			hs ^= hasher(*s);
			++s;
		}
		return hs;
	}

	bool operator()(const char* a, const char* b) const {
		return (std::strcmp(a, b) == 0);
	}
};

// *****************************************************************************************
//
class Murmur64_base {

protected:
	uint64_t load64(const char*& p) const {
		uint64_t x = (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);
		x <<= 8;	x += (uint64_t)(*p++);

		return x;
	}

	uint64_t rotl64(uint64_t x, int8_t r) const 
	{
		return (x << r) | (x >> (64 - r));
	}

	size_t fmix64(uint64_t h) const {
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;
		return h;
	}
};


// *****************************************************************************************
//
template <class string_t>
class Murmur64_simple : public Murmur64_base {
public:

	size_t operator()(const string_t& sv) const {

		size_t h = 0;
		const char* p = sv.data();

		for (int i = 0; i < sv.size() / 8; ++i) {
			uint64_t x = load64(p);
			h ^= fmix64(x);
		}

		uint64_t x = 0;
		switch (sv.size() & 7) { // last 3 bits
		case  7: x |= ((uint64_t)p[6]) << 48; [[fallthrough]];
		case  6: x |= ((uint64_t)p[5]) << 40; [[fallthrough]];
		case  5: x |= ((uint64_t)p[4]) << 32; [[fallthrough]];
		case  4: x |= ((uint64_t)p[3]) << 24; [[fallthrough]];
		case  3: x |= ((uint64_t)p[2]) << 16; [[fallthrough]];
		case  2: x |= ((uint64_t)p[1]) << 8; [[fallthrough]];
		case  1: x |= ((uint64_t)p[0]) << 0;
		}

		h ^= fmix64(x);

		return h;
	}
};


// *****************************************************************************************
//
template <class string_t>
class Murmur64_full : public Murmur64_base {
public:
	size_t operator()(const string_t& sv) const {

		uint64_t h1 = 0;
		uint64_t h2 = 0;

		const uint64_t c1 = 0x87c37b91114253d5ull;
		const uint64_t c2 = 0x4cf5ad432745937full;

		const char* p = sv.data();

		// process 16-byte portions
		int n_portions = (int)(sv.size() / 16);
		for (int i = 0; i < n_portions; ++i) {
			uint64_t k1 = load64(p);
			uint64_t k2 = load64(p);

			k1 *= c1; k1 = rotl64(k1, 31); k1 *= c2; h1 ^= k1;
			h1 = rotl64(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;

			k2 *= c2; k2 = rotl64(k2, 33); k2 *= c1; h2 ^= k2;
			h2 = rotl64(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
		}

		// process the remaining part
		uint64_t k1 = 0;
		uint64_t k2 = 0;

		switch (sv.size() & 15) { // last 4 bits
		case 15: k2 ^= ((uint64_t)p[14]) << 48; [[fallthrough]];
		case 14: k2 ^= ((uint64_t)p[13]) << 40; [[fallthrough]];
		case 13: k2 ^= ((uint64_t)p[12]) << 32; [[fallthrough]];
		case 12: k2 ^= ((uint64_t)p[11]) << 24; [[fallthrough]];
		case 11: k2 ^= ((uint64_t)p[10]) << 16; [[fallthrough]];
		case 10: k2 ^= ((uint64_t)p[9]) << 8; [[fallthrough]];
		case  9: k2 ^= ((uint64_t)p[8]) << 0;
			k2 *= c2; k2 = rotl64(k2, 33); k2 *= c1; h2 ^= k2;
			[[fallthrough]];
		case  8: k1 ^= ((uint64_t)p[7]) << 56; [[fallthrough]];
		case  7: k1 ^= ((uint64_t)p[6]) << 48; [[fallthrough]];
		case  6: k1 ^= ((uint64_t)p[5]) << 40; [[fallthrough]];
		case  5: k1 ^= ((uint64_t)p[4]) << 32; [[fallthrough]];
		case  4: k1 ^= ((uint64_t)p[3]) << 24; [[fallthrough]];
		case  3: k1 ^= ((uint64_t)p[2]) << 16; [[fallthrough]];
		case  2: k1 ^= ((uint64_t)p[1]) << 8; [[fallthrough]];
		case  1: k1 ^= ((uint64_t)p[0]) << 0;
			k1 *= c1; k1 = rotl64(k1, 31); k1 *= c2; h1 ^= k1;
		};

		h1 ^= (uint64_t)sv.size(); h2 ^= (uint64_t)sv.size();

		h1 += h2;
		h2 += h1;

		h1 = fmix64(h1);
		h2 = fmix64(h2);

		h1 += h2;
		h2 += h1;

		return h1 ^ h2;
	}
};
