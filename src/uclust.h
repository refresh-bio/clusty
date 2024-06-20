#pragma once
#include "distances.h"
#include "clustering.h"

#include <vector>
#include <algorithm>

#include <vector>
#include <unordered_map>

template <class DistanceMatrix>
class UClust : public IClustering<DistanceMatrix> {
public:

	int operator()(
		const DistanceMatrix& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override {

		std::unordered_map<int,int> seeds2clusters;
		int n_objects = objects.size();

		assignments.resize(objects.size(), -1);

		// first element always used as a seed
		int obj = objects[0];
		assignments[obj] = 0;
		seeds2clusters[obj] = 0;

		// assign other elements or make them seeds
		for (int i = 1; i < n_objects; ++i) {
			obj = objects[i];

			dist_t max_edge{ 0, dist_t::MAX };
			const dist_t* closest = &max_edge;
			
			// select closest seed among neighbours
			for (const dist_t* edge = distances.begin(obj); edge < distances.end(obj); ++edge) {
				
				if (seeds2clusters.count(edge->u.s.hi) && edge->d < closest->d) {
					closest = edge;
				}
			}
			
			// if distance to seed is below threshold
			if (closest->d <= threshold) {
				assignments[obj] = seeds2clusters[closest->u.s.hi];
			}
			else {
				int cluster_id = seeds2clusters.size();
				seeds2clusters.insert({ obj, cluster_id });
				assignments[obj] = cluster_id;
			}
		}

		return seeds2clusters.size();
	}
};


