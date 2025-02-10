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
class UClust : public IClustering<Distance> {
public:

	int operator()(
		SparseMatrix<Distance>& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override {

		std::unordered_map<int,int> seeds2clusters;
		int n_objects = (int)objects.size();

		assignments.resize(objects.size(), -1);

		// first element always used as a seed
		int obj = objects[0];
		assignments[obj] = 0;
		seeds2clusters[obj] = 0;

		// assign other elements or make them seeds
		for (int i = 1; i < n_objects; ++i) {
			obj = objects[i];

			Distance max_edge;
			const Distance* closest = &max_edge;
			
			// select closest seed among neighbours
			for (const Distance* edge = distances.begin(obj); edge < distances.end(obj); ++edge) {
				
				if (seeds2clusters.count(edge->get_id()) && edge->get_d() < closest->get_d()) {
					closest = edge;
				}
			}
			
			// if distance to seed is below threshold
			if (closest->get_d() <= threshold) {
				assignments[obj] = seeds2clusters[closest->get_id()];
			}
			else {
				int cluster_id = (int)seeds2clusters.size();
				seeds2clusters.insert({ obj, cluster_id });
				assignments[obj] = cluster_id;
			}
		}

		return (int)seeds2clusters.size();
	}
};


