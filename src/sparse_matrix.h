// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "distances.h"

#include <vector>


/*********************************************************************************************************************/
template <class Distance>
class SparseMatrix : public IMatrix
{
public:

	size_t n_elements{ 0 };

	std::vector<std::vector<Distance>> distances;

	SparseMatrix() {}

	virtual ~SparseMatrix() {}

	size_t num_objects() const { return distances.size(); }

	size_t num_elements() const { return n_elements; }

	size_t num_neighbours(int i) const { return distances[i].size(); }

	const Distance* begin(int row_id) const { return distances[row_id].data(); }
	const Distance* end(int row_id) const { return distances[row_id].data() + distances[row_id].size(); }

	void clear_row(int row_id) { std::vector<Distance>().swap(distances[row_id]); }
};

