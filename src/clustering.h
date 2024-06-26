#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <vector>
#include <numeric>

#include "distances.h"

struct node_t {
	int first = -1;
	int second = -1;
	double distance = -1.0;

	node_t() {}
	node_t(int first, int second, double distance) 
		: first(first), second(second), distance(distance) {}
};

template <class DistanceMatrix>
class IClustering {
public:

	virtual int operator()(
		const DistanceMatrix& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) = 0;

	virtual ~IClustering() {};

};

template <class DistanceMatrix>
class HierarchicalClustering : public IClustering<DistanceMatrix> {
protected:
	
	void makeDendrogram(
		const std::vector<dist_t>& lambda,
		const std::vector<int>& pi,
		std::vector<node_t>& dendrogram) 
	{
		int n_objects = (int)lambda.size();
		
		std::vector<int> elements(n_objects - 1);
		std::iota(elements.begin(), elements.end(), 0);

		stable_sort(elements.begin(), elements.end(), [&lambda](int x, int y) {
			return lambda[x] < lambda[y];
			});

		std::vector<int> index(n_objects);
		std::iota(index.begin(), index.end(), 0);

		dendrogram.reserve(2 * n_objects);
		dendrogram.resize(n_objects);
		
		for (int i = 0; i < n_objects - 1; ++i) {
			int j = elements[i];
			int next = pi[j];
			dendrogram.emplace_back(index[j], index[next], lambda[j].d);
			index[next] = n_objects + i;
		}
	}
	
	int dendrogramToAssignments(
		const std::vector<node_t>& dendrogram,
		double threshold,
		std::vector<int>& assignments) {

		int n_objects = assignments.size();

		std::vector<int> prevs(dendrogram.size() + 1, -1);
		std::vector<int> num_visits(dendrogram.size() + 1, 0);
		
		int cluster_id = 0;
		
		// dendrogram can have multiple isolated roots - iterate over positions from the last one
		for (int last_pos = dendrogram.size() - 1; last_pos >= 0; --last_pos) {

			// visited nodes are inner roots - omit them
			if (num_visits[last_pos] > 0) {
				continue;
			}

			int cur_pos = last_pos;
			while (true) {

				if (cur_pos < n_objects) {
					// sequence was reached
					++num_visits[cur_pos];
					assignments[cur_pos] = cluster_id;

					// singleton sequence
					if (cur_pos == last_pos) {
						break;
					}
					
					cur_pos = prevs[cur_pos];
				}
				else {
					// internal node

					if (num_visits[cur_pos] == 0) {
						// no visits - left branch
						int dest_pos = dendrogram[cur_pos].first;
						++num_visits[cur_pos];
						prevs[dest_pos] = cur_pos;
						cur_pos = dest_pos;

					}
					else if (num_visits[cur_pos] == 1) {
						// one visit - right branch

						// if distance between nodes is larger than the threshold - assign to different cluster
						if (dendrogram[cur_pos].distance > threshold) {
							++cluster_id;
						}

						int dest_pos = dendrogram[cur_pos].second;
						++num_visits[cur_pos];
						prevs[dest_pos] = cur_pos;
						cur_pos = dest_pos;
					}
					else {
						// two visits - node processed
						if (cur_pos == last_pos) {
							// root processed

							break;
						}

						++num_visits[cur_pos];
						cur_pos = prevs[cur_pos];
					}
				}

			}

			++cluster_id; // next potential root - next cluster
		}

		return cluster_id;
	}
};
