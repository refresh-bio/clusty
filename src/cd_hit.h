#pragma once
#include "distances.h"
#include "clustering.h"

#include <vector>
#include <algorithm>

#include <vector>
#include <unordered_map>

template <class DistanceMatrix>
class CdHit : public IClustering<DistanceMatrix> {
public:

	int operator()(
		const DistanceMatrix& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override {

		int n_objects = objects.size();

		assignments.resize(objects.size(), -1);

		int cluster_id = 0;

		for (int i = 0; i < n_objects; ++i) {
			auto obj = objects[i];

			if (assignments[obj] == -1) {
				// unassigned object - make a new seed
				assignments[obj] = cluster_id;

				// iterate over connected object and assign those which are unassigned
				for (const dist_t* edge = distances.begin(obj); edge < distances.end(obj); ++edge) {
					int other = edge->u.s.hi;

					if (edge->d <= threshold && assignments[other] == -1) {
						assignments[other] = cluster_id;
					}	
				}

				++cluster_id;
			}
		}

		return cluster_id;
	}
};


