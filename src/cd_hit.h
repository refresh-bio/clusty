#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "distances.h"
#include "clustering.h"

#include <vector>
#include <algorithm>

#include <vector>
#include <unordered_map>

template <class Distance>
class CdHit : public IClustering<Distance> {
public:

	int operator()(
		SparseMatrix<Distance>& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override {

		int n_objects = (int)objects.size();

		assignments.resize(objects.size(), -1);

		int cluster_id = 0;

		for (int i = 0; i < n_objects; ++i) {
			auto obj = objects[i];

			if (assignments[obj] == -1) {
				// unassigned object - make a new seed
				assignments[obj] = cluster_id;

				// iterate over connected object and assign those which are unassigned
				for (const Distance* edge = distances.begin(obj); edge < distances.end(obj); ++edge) {
					int other = edge->get_id();

					if (edge->get_d() <= threshold && assignments[other] == -1) {
						assignments[other] = cluster_id;
					}	
				}

				++cluster_id;
			}
		}

		return cluster_id;
	}
};
