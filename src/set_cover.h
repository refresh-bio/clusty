#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <vector>
#include <algorithm>
 
#include "distances.h"
#include "clustering.h"

#define debug(x) std::cerr << __FILE__ << " (" << __LINE__ << ") " << #x << " == " << (x) << std::endl


/** MMseqs0 clustering algorithm */
template <class Distance>
class SetCover : public IClustering<Distance> 
{
public:
	int operator()(
		SparseMatrix<Distance>& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override
	{
		
		const int NO_ASSIGNMENT { -1 };
		int nObjects = (int)objects.size();
		assignments = std::vector<int> (nObjects, NO_ASSIGNMENT); // -1: no assignment
		
		std::vector<std::pair<int, int>> obj2connections(nObjects);
		for (int i = 0; i < nObjects; ++i) {
			int obj = objects[i];
			obj2connections[i].first = obj;
			obj2connections[i].second = distances.num_neighbours(obj);
		}

		std::stable_sort(obj2connections.begin(), obj2connections.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
		
 		int cluster_number = 0;
		// We analyse nodes from the one with the greatest number of neighbours.
 		for (int i = 0; i < nObjects; ++i)
 		{
			int obj = obj2connections[i].first;
			
 			if (assignments[obj] == NO_ASSIGNMENT)
 			{
				assignments[obj] = cluster_number; // seed of a new cluster ... 
				
				// ... and its neighbours:
				for (const Distance* edge = distances.begin(obj); edge < distances.end(obj); ++edge) {
					auto other = edge->get_id();
					if (edge->get_d() <= threshold && assignments[other] == NO_ASSIGNMENT) {
						assignments[other] = cluster_number;
					}
				}
 				cluster_number++;
 			}
 		}
		return cluster_number;		
	}
	
};
