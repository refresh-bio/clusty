// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "graph_sparse.h"
#include "hasher.h"
#include "io.h"
#include "chunked_vector.h"

#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <string>

// *******************************************************************************************/
union name_or_id_t {
	std::string_view name;
	int id;

	name_or_id_t() : id{0} {}
};

using NamedEdgesCollection = EdgesCollection<name_or_id_t>;

/*********************************************************************************************************************/
template <class Distance>
class GraphNamed : public GraphSparse<Distance> {

	using ids_pair_t = std::pair<int, int>;

	std::unordered_map<std::string_view, ids_pair_t, Murmur64_full<std::string_view>> names2ids;

	std::vector<std::string_view> ids2names;

	chunked_vector<char> namesBuffer{ 16LL << 20 }; // 16MB chunk size

public:
	GraphNamed(int numThreads) : GraphSparse<Distance>(numThreads) {}

	size_t getNumInputVertices() const override { return this->names2ids.size(); }

	void reorderObjects(
		const std::vector<std::string_view>& externalNames,
		std::vector<int>& objects) const override {

		int obj_id = 0;
		for (const auto& name : externalNames) {
			int local_id = get_id(name);
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

	void print(std::ostream& out) const override;

protected:

	int get_id(const std::string_view& name) const {
		auto it = names2ids.find(name);
		if (it == names2ids.end()) {
			return -1;
		}
		else {
			return it->second.first;
		}
	}

	IEdgesCollection* createEdgesCollection(size_t preallocSize) override {
		return new NamedEdgesCollection(preallocSize);
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
void GraphNamed<Distance>::initLoad() {

	// invoke superclass method
	GraphSparse<Distance>::initLoad();

	// assume space for 8M objects
	ids2names.reserve(8LL << 20);

	name_or_id_t ai;
}


/*********************************************************************************************************************/
template <class Distance>
bool GraphNamed<Distance>::parseBlock(
	char* block_begin,
	char* block_end,
	distance_transformation_t transform,
	IEdgesCollection& edges,
	size_t& n_rows) {

	NamedEdgesCollection& namedEdges = dynamic_cast<NamedEdgesCollection&>(edges);

	int n_columns = (int)this->filters.size();
	n_rows = 0;

	char* line = block_begin;

	while (line != block_end) {

		++n_rows;

		char* p = line;
		NamedEdgesCollection::edge_t edge;
		edge.second = std::numeric_limits<double>::max();
		bool carryOn = true;
		int k = 0;
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

			size_t name_len = q - p;

			if (k < 2 && c == this->sequenceColumnIds[k]) {
				edge.first[k].name = std::string_view(p, name_len);
				++k;
			}
			else if (c == this->distanceColumnId || this->filters[c].enabled) {
				double value = Conversions::strtod(p, &p);

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

		// Push edge if it passed filters (self-loops will be filtered in updateMatrix,
		// but we need to push them here to register the node names)
		if (carryOn) {
			namedEdges.data.push_back(edge);
		}

		// if new line was detected
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

	return false; // cannot free input buffer
}


/*********************************************************************************************************************/
template <class Distance>
void GraphNamed<Distance>::updateMappings(
	IEdgesCollection& edges) {

	NamedEdgesCollection& namedEdges{ dynamic_cast<NamedEdgesCollection&>(edges) };

	for (NamedEdgesCollection::edge_t& e : namedEdges.data) {

		decltype(names2ids.begin()) its[2];

		for (int k = 0; k < 2; ++k) {
			const std::string_view& name = e.first[k].name;

			// store name in hashtable
			its[k] = names2ids.find(name);

			if (its[k] == names2ids.end()) {

				char* dst = namesBuffer.resize_for_additional(name.size() + 1);
				std::copy_n(name.data(), name.size(), dst); // 0 is already there
				auto it_and_flag = names2ids.insert({ std::string_view(dst, name.size()), {-1, names2ids.size()} }); // -1 indicate singleton

				its[k] = it_and_flag.first;
			}

			auto it = its[k];

			// if name not mapped to numerical ids
			if (it->second.first == -1) {
				ids2names.push_back(it->first);
				it->second.first = (int)ids2names.size() - 1;
			}

			e.first[k].id = it->second.first;
		}
	}
}

/*********************************************************************************************************************/
template <class Distance>
void GraphNamed<Distance>::extendMatrix() {
	// can be only larger
	this->matrix.distances.resize(ids2names.size());
}

/*********************************************************************************************************************/
template <class Distance>
void GraphNamed<Distance>::updateMatrix(
	const IEdgesCollection& edges,
	int startRow,
	int stride) {

	const NamedEdgesCollection& namedEdges{ dynamic_cast<const NamedEdgesCollection&>(edges) };

	for (const NamedEdgesCollection::edge_t& e : namedEdges.data) {
		// Skip self-loops (diagonal elements)
		if (e.first[0].id == e.first[1].id) {
			continue;
		}

		for (int k = 0; k < 2; ++k) {

			int lid = e.first[k].id;
			if ((lid % stride) == startRow && e.second < std::numeric_limits<double>::max()) {
				auto& D = this->matrix.distances[lid];

				// extend capacity by factor 1.5 with 16 as an initial state
				if (D.capacity() == D.size()) {
					D.reserve(D.capacity() == 0 ? 16 : size_t(D.capacity() * 1.5));
				}

				D.emplace_back(e.first[k ^ 1].id, e.second);
			}
		}
	}
}


/*********************************************************************************************************************/
template <class Distance>
int GraphNamed<Distance>::saveAssignments(
	std::ofstream& ofs,
	const std::vector<std::string_view>& globalNames,
	const std::vector<int>& assignments,
	char separator,
	bool useRepresentatives) const {

	std::vector<int> old2new;
	this->sortClustersBySize(assignments, old2new);

	int singleton_id = (int)old2new.size();
	if (globalNames.empty()) {

		std::vector<std::tuple<std::string_view, int>> names_n_clusters(assignments.size());

		int i = 0;
		std::transform(assignments.begin(), assignments.end(), names_n_clusters.begin(), [this, &i, &old2new](int a) {
			return std::make_tuple(this->ids2names[i++], old2new[a]);
			});

		std::sort(names_n_clusters.begin(), names_n_clusters.end(), [](const auto& p, const auto& q) {
			return (std::get<1>(p) == std::get<1>(q)) ? (std::get<0>(p) < std::get<0>(q)) : (std::get<1>(p) < std::get<1>(q));
			});


		if (useRepresentatives) {
			std::vector<std::tuple<std::string_view, std::string_view>> names_n_reps;
			this->fillRepresentatives(names_n_clusters, names_n_reps);
			saveTableBuffered<2>(ofs, std::array<std::string, 2>({ "object", "cluster" }), names_n_reps, separator);

		}
		else {
			saveTableBuffered<2>(ofs, std::array<std::string, 2>({ "object", "cluster" }), names_n_clusters, separator);
		}

	}
	else {

		std::vector<std::tuple<std::string_view, int, int>> names_n_clusters_n_ids(globalNames.size());

		int inside_id = 0;
		int outsize_id = (int)assignments.size();

		for (int gi = 0; gi < globalNames.size(); ++gi) {
			auto& name = globalNames[gi];

			int local_id = get_id(name);
			if (local_id == -1) {
				// not in matrix
				if (outsize_id >= names_n_clusters_n_ids.size()) {
					throw std::runtime_error("Names mismatch between distance and objects files.");
				}

				auto& out{ names_n_clusters_n_ids[outsize_id++] };
				std::get<0>(out) = name;
				std::get<1>(out) = singleton_id++;
				std::get<2>(out) = gi;
			}
			else {
				auto& out{ names_n_clusters_n_ids[inside_id++] };
				// in matrix
				std::get<0>(out) = name;
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

	return singleton_id;
}

/*********************************************************************************************************************/
template <class Distance>
void GraphNamed<Distance>::print(std::ostream& out) const {

	out.precision(10);

	std::vector<std::string_view> names(names2ids.size());
	int i = 0;
	for (auto q : names2ids) {
		names[i] = q.first;
		++i;
	}

	std::sort(names.begin(), names.end());

	for (auto name : names) {

		int i = names2ids.at(name).first;
		const std::vector<Distance>& row = this->matrix.distances[i];

		for (auto& p : row) {
			out << ids2names[i] << "," << ids2names[p.get_id()] << "," << std::fixed << p.get_d() << std::endl;
		}
	}
}

