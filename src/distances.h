#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#include "conversion.h"

#include <limits>
#include <algorithm>
#include <fstream>
#include <cstdint>
#include <functional>


/*********************************************************************************************************************/
using distance_transformation_t = std::function<double(double)>;

/*********************************************************************************************************************/
#pragma pack(push, 4)
class dist_t {
	
	double d;		///< distance   std::numeric_limits<double>::max() 
	uint32_t id;	// object identifier

public:
	dist_t() : d(std::numeric_limits<double>::max()), id(0) {}
	dist_t(uint32_t id, double d) : d(d), id(id) {}

	double get_d() const { return d; }
	uint32_t get_id() const { return id; }

	bool operator<(const dist_t& rhs) const {
		return (id == rhs.id) ? (d < rhs.d) : (id < rhs.id);
	}
};
#pragma pack(pop)

/*********************************************************************************************************************/
class mini_dist_t {
	uint32_t id;	// object identifier

public:
	mini_dist_t() : id(0) {}
	mini_dist_t(uint32_t id, double d) : id(id) {}

	double get_d() const { return 0.0; } // mini_dist_t is always below threshold
	uint32_t get_id() const { return id; }

	bool operator<(const mini_dist_t& rhs) const {
		return (id < rhs.id);
	}
};

/*********************************************************************************************************************/
class IMatrix {
public:
	virtual ~IMatrix() {}
};


