// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "graph_sparse.h"
#include "log.h"
#include "io.h"
#include "chunked_vector.h"
#include "conversion.h"

#include <vector>
#include <limits>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <array>


// *******************************************************************************************/
using NumberedEdgesCollection = EdgesCollection<int>;

// *******************************************************************************************/
template <class Distance>
class GraphNumbered : public GraphSparse<Distance> {

	std::vector<int> global2local;
	std::vector<int> local2global;

public:

	GraphNumbered(int numThreads) : GraphSparse<Distance>(numThreads) {}

	size_t getNumInputVertices() const override { return local2global.size(); }

	void reorderObjects(
		const std::vector<std::string_view>& externalNames,
		std::vector<int>& objects) const override {
	
		int obj_id = 0;
		for (int i = 0; i < (int)externalNames.size(); ++i) {
			int local_id = get_local_id(i);
			if (local_id != -1) {
				objects[obj_id++] = local_id;
			}
		}
	}

	int saveAssignments(
		std::ofstream& ofs,
		const std::vector<std::string_view>& globalNames,
		const std::vector<int>& assignments,
		char separator,
		bool useRepresentatives) const override;

	void print(std::ostream& out) const override {}



protected:

	int get_local_id(int global_id) const {
		if ((size_t)global_id > global2local.size() - 1) {
			return -1;
		}
		else {
			return global2local[global_id];
		}
	}

	IEdgesCollection* createEdgesCollection(size_t preallocSize) override {
		return new NumberedEdgesCollection(preallocSize);
	};

	void initLoad() override;

	bool parseBlock(
		char* block_begin,
		char* block_end,
		distance_transformation_t transform,
		IEdgesCollection& edges,
		size_t& n_rows)  override;

	void updateMappings(
		IEdgesCollection& edges) override;

	void extendMatrix() override;

	void updateMatrix(
		const IEdgesCollection& edges,
		int startRow,
		int stride) override;

};


/*********************************************************************************************************************/
template <class Distance>
void GraphNumbered<Distance>::initLoad() {

	// invoke superclass method
	GraphSparse<Distance>::initLoad();
	
	// assume space for 8M objects
	global2local.reserve(8LL << 20);
	local2global.reserve(8LL << 20);
}

/*********************************************************************************************************************/
template <class Distance>
bool GraphNumbered<Distance>::parseBlock(
	char* block_begin,
	char* block_end,
	distance_transformation_t transform,
	IEdgesCollection& edges,
	size_t& n_rows) {
	
	auto& numEdgesCol{ dynamic_cast<NumberedEdgesCollection&>(edges) };
	auto& v_edges { numEdgesCol.data };
	
	int n_columns = (int)this->filters.size();
	n_rows = 0;

	char* line = block_begin;
	
	int max_value = 0;

	while (line != block_end) {
		
		++n_rows;

		char* p = line;
		int k = 0;
		
		NumberedEdgesCollection::edge_t edge;
		edge.second = std::numeric_limits<double>::max();
		bool carryOn = true;
		bool reachedNewline = false;

		for (int c = 0; c < n_columns; ++c) {
			char* q = std::find_if(p, block_end, this->isSeparator); // support both tsv and csv files
		
			reachedNewline = this->isNewline(*q);

			if (q != block_end) {
				if ((c == n_columns - 1) ^ reachedNewline) {
					throw std::runtime_error("Ill-formatted input table in row " + std::to_string(n_rows) + "\n" + std::string(line, line + 50));
				}
				*q = 0;
			}

			if (k < 2 && c == this->sequenceColumnIds[k]) {
				edge.first[k] = Conversions::strtol(p, nullptr);
				++k;
			}
			else if (c == this->distanceColumnId || this->filters[c].enabled) {
				double value = Conversions::strtod(p, nullptr);

				if (c == this->distanceColumnId) {
					edge.second = transform(value); 	// convert similarity to distance if neccessary
				}

				// check distance condition
				if (value < this->filters[c].min || value > this->filters[c].max) {
					p = q + 1;
					//edge.second = std::numeric_limits<double>::max();
					carryOn = false;
					break;
				}
			}

			p = q + 1;
		}

		// do not consider rows diagonal elements (they are assumed to have 0 distance)
		if (carryOn && (edge.first[0] != edge.first[1])) {

			v_edges.push_back(edge);

			if (edge.first[0] > max_value) {
				max_value = edge.first[0];
			}

			if (edge.first[1] > max_value) {
				max_value = edge.first[1];			
			}
		}

		// if new line was detected or end of block reached
		if (reachedNewline || p > block_end) {
			// decrease p so it points 0 which replaced the newline
			--p;	
		}
		else {
			// find new line character
			p = std::find_if(p, block_end, this->isNewline);
		}
		
		// p should be at symbol right after newline (but new line can consists of two chars)
		line = std::find_if(p, block_end, [](char c) { return c != '\r' && c != '\n' && c != 0; });
	}

	// make the largest id at index 1 of the front element 
	numEdgesCol.maxEdge = max_value;

	return true;
}


/*********************************************************************************************************************/
template <class Distance>
void GraphNumbered<Distance>::updateMappings(
	IEdgesCollection& edges) {

	NumberedEdgesCollection& numberedEdges{ dynamic_cast<NumberedEdgesCollection&>(edges) };

	// adjust size if needed
	if (numberedEdges.maxEdge >= (int)global2local.size()) {
		global2local.resize(numberedEdges.maxEdge + 1, -1);
	}

	// update mappings in collection
	for (NumberedEdgesCollection::edge_t& e : numberedEdges.data) {
		
		for (int k = 0; k < 2; ++k) {

			int gid = e.first[k];
			int lid = global2local[gid];

			if (lid == -1) {
				lid = (int)local2global.size();
				local2global.push_back(gid);
				global2local[gid] = lid;
			}

			e.first[k] = lid;
		}
	}

	
}


/*********************************************************************************************************************/
template <class Distance>
void GraphNumbered<Distance>::extendMatrix() {
	// can be only larger
	this->matrix.distances.resize(local2global.size());
}


/*********************************************************************************************************************/
template <class Distance>
void GraphNumbered<Distance>::updateMatrix(
	const IEdgesCollection& edges,
	int startRow,
	int stride) {

	const NumberedEdgesCollection& numberedEdges{ dynamic_cast<const NumberedEdgesCollection&>(edges) };

	for (const NumberedEdgesCollection::edge_t& e : numberedEdges.data) {

		for (int k = 0; k < 2; ++k) {
			int lid = e.first[k];

			if ((lid % stride) == startRow && e.second < std::numeric_limits<double>::max()) {
				auto& D = this->matrix.distances[lid];
				
				// extend capacity by factor 1.5 with 16 as an initial state
				if (D.capacity() == D.size()) {
					D.reserve(D.capacity() == 0 ? 16 : size_t(D.capacity() * 1.5));
				}

				D.emplace_back(e.first[k ^ 1], e.second);
			}
		}
	}
}


/*********************************************************************************************************************/
template <class Distance>
int GraphNumbered<Distance>::saveAssignments(
	std::ofstream& ofs,
	const std::vector<std::string_view>& globalNames,
	const std::vector<int>& assignments,
	char separator,
	bool useRepresentatives) const {

	std::vector<int> old2new;
	this->sortClustersBySize(assignments, old2new);

	int singleton_id = (int)old2new.size();

	if (globalNames.empty()) {
		std::vector<std::tuple<int, int>> ids_n_clusters(assignments.size());

		int i = 0;
		std::transform(assignments.begin(), assignments.end(), ids_n_clusters.begin(), [this, &i, &old2new](int a) {
			return std::make_tuple(this->local2global[i++], old2new[a]);
		});

		// sort increasingly by cluster and by object id 
		std::sort(ids_n_clusters.begin(), ids_n_clusters.end(), [](const auto& p, const auto& q) {
			return (std::get<1>(p) == std::get<1>(q)) ? (std::get<0>(p) < std::get<0>(q)) : (std::get<1>(p) < std::get<1>(q));
		});

		if (useRepresentatives) {
			std::vector<std::tuple<int, int>> ids_n_reps;
			this->fillRepresentatives(ids_n_clusters, ids_n_reps);
			saveTableBuffered<2>(ofs, std::array<std::string, 2>({ "object", "cluster" }), ids_n_reps, separator);
		}
		else {
			saveTableBuffered<2>(ofs, std::array<std::string, 2>({ "object", "cluster" }), ids_n_clusters, separator);
		}
	}
	else {
		
		std::vector<std::tuple<std::string_view, int, int>> names_n_clusters_n_ids(globalNames.size());

		int inside_id = 0;
		int outsize_id = (int)assignments.size();
		
		for (int gi = 0; gi < globalNames.size(); ++gi) {
		
			int local_id = get_local_id(gi);
			if (local_id == -1) {
				// not in matrix 
				auto& out{ names_n_clusters_n_ids[outsize_id++] };
				std::get<0>(out) = globalNames[gi];
				std::get<1>(out) = singleton_id++;
				std::get<2>(out) = gi;
			}
			else {
				auto& out{ names_n_clusters_n_ids[inside_id++] };
				// in matrix
				std::get<0>(out) = globalNames[gi];
				std::get<1>(out) = old2new[assignments[local_id]];
				std::get<2>(out) = gi;
			}
		}

		// sort increasingly inside part by cluster and by object id 
		std::sort(names_n_clusters_n_ids.begin(), names_n_clusters_n_ids.begin() + inside_id, [](const auto& p, const auto& q) {
			return (std::get<1>(p) == std::get<1>(q)) ? (std::get<2>(p) < std::get<2>(q)) : (std::get<1>(p) < std::get<1>(q));
		});

		if (useRepresentatives) {
			std::vector<std::tuple<std::string_view, std::string_view>> names_n_reps;
			this->fillRepresentatives(names_n_clusters_n_ids, names_n_reps);
			saveTableBuffered<2>(ofs, std::array<std::string, 2>({ "object", "cluster" }), names_n_reps, separator);

		}
		else {
			saveTableBuffered<2>(ofs, std::array<std::string, 2>({ "object", "cluster" }), names_n_clusters_n_ids, separator);
		}
	}

	return singleton_id; // return total number of clusters
}

