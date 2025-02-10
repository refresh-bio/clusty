#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "clustering.h"
#include <stdexcept>


#ifndef NO_LEIDEN
#include <igraph.h>
#endif

struct LeidenParams {

	double resolution{ 0.7 };
	double beta{ 0.01 };
	int numIterations{ 2 };
};


template <class Distance>
class Leiden : public IClustering<Distance> {

private:
	LeidenParams params;

public:


#ifdef NO_LEIDEN

	Leiden(const LeidenParams& params) : params(params) {
		throw std::runtime_error("Leiden algorithm not available.");
	}

	int operator()(
		SparseMatrix<Distance>& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override {

		throw std::runtime_error("Leiden algorithm not available.");
		return 0;
	}
#else
	
	Leiden(const LeidenParams& params) : params(params) {}
	
	int operator()(
		SparseMatrix<Distance>& distances,
		const std::vector<int>& objects,
		double threshold,
		std::vector<int>& assignments) override {


		igraph_t g;
		igraph_vector_t edge_weights;

		load_graph(distances, g, edge_weights);

		igraph_integer_t n_clusters;
		igraph_vector_int_t memb_vec;

		igraph_vector_int_init(&memb_vec, distances.num_objects());
		igraph_community_leiden(&g, &edge_weights, nullptr,
			params.resolution, params.beta, false, params.numIterations,
			&memb_vec, &n_clusters, nullptr);

		assignments.resize(distances.num_objects(), -1);

		for (size_t i = 0; i < distances.num_objects(); ++i) {
			assignments[i] = igraph_vector_int_get(&memb_vec, i);
		}

		return (int)n_clusters;
	}


	void load_graph(SparseMatrix<Distance>& matrix, igraph_t& g, igraph_vector_t& edge_weights) {

		igraph_vector_int_t edges;
		igraph_vector_int_init(&edges, 0);

		igraph_vector_init(&edge_weights, 0);

		for (int i = 0; i < matrix.num_objects(); ++i) {
			for (const Distance* edge = matrix.begin(i); edge < matrix.end(i); ++edge) {
				igraph_vector_int_push_back(&edges, i);
				igraph_vector_int_push_back(&edges, edge->get_id());
				igraph_vector_push_back(&edge_weights, 1.0 - edge->get_d());
			}

			matrix.clear_row(i);
		}

		/*
		for (const dist_t& edge : matrix.get_distances()) {
			igraph_vector_int_push_back(&edges, edge.u.s.lo);
			igraph_vector_int_push_back(&edges, edge.u.s.hi);
			igraph_vector_push_back(&edge_weights, 1.0 - edge.d); 
			//igraph_vector_push_back(&edge_weights, 1);
		}
		*/

		igraph_empty(&g, matrix.num_objects(), IGRAPH_UNDIRECTED);
		igraph_add_edges(&g, &edges, nullptr);
	}
#endif

};