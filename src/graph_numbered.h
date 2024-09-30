// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "sparse_matrix.h"
#include "graph.h"
#include "log.h"

#include <vector>
#include <limits>
#include <unordered_map>
#include <iostream>
#include <fstream>


template <class Distance>
class GraphNumbered : public Graph {

	std::vector<int> global2local;
	std::vector<int> local2global;

	SparseMatrix<Distance> matrix;

public:

	IMatrix& getMatrix() override { return matrix; }

	size_t getNumVertices() const override { return matrix.num_objects(); }

	size_t getNumInputVertices() const override { return global2local.size(); }

	size_t getNumEdges() const override { return matrix.num_elements(); }

	void reorderObjects(
		const std::vector<std::string>& externalNames,
		std::vector<int>& objects) const override {
	
		int obj_id = 0;
		for (size_t i = 0; i < externalNames.size(); ++i) {
			int local_id = get_local_id(i);
			if (local_id != -1) {
				objects[obj_id++] = local_id;
			}
		}
	}

	size_t load(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		distance_transformation_t transform,
		const std::map<std::string, ColumnFilter>& columns2filters) override;

	int saveAssignments(
		std::ofstream& ofs,
		const std::vector<std::string>& globalNames,
		const std::vector<int>& assignments,
		char separator) const override;

	int saveRepresentatives(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator) const override;

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
};



/*********************************************************************************************************************/
template <class Distance>
size_t GraphNumbered<Distance>::load(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	distance_transformation_t transform,
	const std::map<std::string, ColumnFilter>& columns2filters) {

	size_t n_total_distances = 0;
	const size_t N = 512ULL << 20;
	char* buf = new char[N];
	char* buf_end = buf + N - 1;

	// get header
	int col_ids[]{ 0, 1 };
	int col_distance = 2; // by default use 3rd column as the one with distance 
	std::vector<ColumnFilter> filters;
	processHeader(ifs, idColumns, distanceColumn, columns2filters, col_ids, col_distance, filters);
	int n_columns = filters.size();

	auto is_sep = [](char c) {return c == ',' || c == '\t' || c == '\r' || c == '\t'; };
	auto is_newline = [](char c) {return c == '\r' || c == '\n'; };

	// assume space for 8M objects
	matrix.distances.reserve(8LL << 20);
	global2local.reserve(8LL << 20);
	local2global.reserve(8LL << 20);

	bool continueReading = true;
	char* place = buf;

	while (continueReading) {
		size_t n_wanted = buf_end - place;
		ifs.read(place, n_wanted);
		size_t n_read = ifs.gcount();
		char* block_end = place + n_read;
		int offset = 0;

		// no more data
		if (n_read < n_wanted) {
			continueReading = false;
		}
		else {
			// find last newline
			while (!is_newline(*(block_end - 1))) {
				--block_end;
				++offset;
			}
		}

		//LOG_DEBUG << "portion: " << n_read << ", offset: " << offset << ", carryOn: " << continueReading << endl;

		// pass through buffer
		char* line = buf;
		while (true) {
			char* line_end = std::find_if(line, block_end, is_newline);

			// no more lines
			if (line_end == block_end) {
				break;
			}

			++n_total_distances;

			char* p = line;
			double d = std::numeric_limits<double>::max();
			int k = 0;
			bool carryOn = true;

			int global_ids[2];

			for (int c = 0; c < n_columns; ++c) {
				char* q = std::find_if(p, line_end, is_sep); // support both tsv and csv files
				*q = 0;

				if (k < 2 && c == col_ids[k]) {
					global_ids[k] = Conversions::strtol(p, nullptr);

					if (global_ids[k] + 1 > (int)global2local.size()) {
						global2local.resize(global_ids[k] + 1, -1);
					}

					++k;
				}
				else if (c == col_distance || filters[c].enabled) {
					double value = Conversions::strtod(p, &p);

					if (c == col_distance) {
						d = transform(value); 	// convert similarity to distance if neccessary
					}

					// check distance condition
					if (value < filters[c].min || value > filters[c].max) {
						carryOn = false;
						break;
					}
				}

				p = q + 1;
			}

			// move to the next line
			line = std::find_if(line_end, block_end, [](char c) { return c != '\r' && c != '\n' && c != 0; });

			// do not consider rows not fulfilling conditions
			if (carryOn == false) {
				continue;
			}

			int pair_ids[2];


			for (int k = 0; k < 2; ++k) {

				int gid = global_ids[k];
				int lid = global2local[gid];

				if (lid == -1) {
					lid = local2global.size();
					local2global.push_back(gid);
					global2local[gid] = lid;
				}

				pair_ids[k] = lid;
			}

			uint32_t i = std::min(pair_ids[0], pair_ids[1]);
			uint32_t j = std::max(pair_ids[0], pair_ids[1]);

			// omit diagonal elements - they are assumed to have 0 distance
			if (i == j) {
				continue;
			}

			if (matrix.distances.size() <= j) {
				matrix.distances.resize(j + 1);
			}

			auto& Di = matrix.distances[i];
			auto& Dj = matrix.distances[j];

			// extend capacity by factor 1.5 with 16 as an initial state
			if (Di.capacity() == Di.size()) {
				Di.reserve(Di.capacity() == 0 ? 16 : size_t(Di.capacity() * 1.5));
			}

			if (Dj.capacity() == Dj.size()) {
				Dj.reserve(Dj.capacity() == 0 ? 16 : size_t(Dj.capacity() * 1.5));
			}

			Di.emplace_back(j, d);
			Dj.emplace_back(i, d);
		}

		// copy remaining part after consuming all the lines
		if (continueReading && offset > 0) {
			memcpy(buf, block_end, offset);
			place = buf + offset;
		}
		else {
			place = buf;
		}
	}

	// if neccessary, sort distances in rows according to the second id
	matrix.n_elements = 0;

	for (auto& row : matrix.distances) {
		std::sort(row.begin(), row.end());
		auto newEnd = std::unique(row.begin(), row.end(), [](const Distance& a, const Distance& b) { return a.get_id() == b.get_id(); });

		row.erase(newEnd, row.end());
		matrix.n_elements += row.size();
	}

	delete[] buf;

	// Print distance histogram in the verbose mode
	if (Log::getInstance(Log::LEVEL_VERBOSE).isEnabled()) {

		std::vector<double> histo_bounds{ 0 };
		double width = 0.001;

		while (histo_bounds.back() < 0.05)
		{
			histo_bounds.push_back(histo_bounds.back() + width);
		}
		histo_bounds.push_back(std::numeric_limits<double>::max());
		std::vector<int> histo(histo_bounds.size());

		for (auto& row : matrix.distances) {
			for (const auto& e : row) {
				for (size_t i = 0; i < histo_bounds.size(); ++i) {
					if (e.get_d() < histo_bounds[i]) {
						++histo[i];
						break;
					}
				}
			}
		}

		LOG_VERBOSE << std::endl << "Distance histogram" << std::endl;
		for (size_t i = 0; i < histo_bounds.size(); ++i) {
			LOG_VERBOSE << "  d < " << histo_bounds[i] << ": " << histo[i] << std::endl;
		}
		LOG_VERBOSE << std::endl;
	}

	return n_total_distances;
}


/*********************************************************************************************************************/
template <class Distance>
int GraphNumbered<Distance>::saveRepresentatives(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator) const
{

	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;

	std::vector<std::string> names;
	if (globalNames.empty()) {
		names.resize(global2local.size());
		int i = 0;
		std::generate(names.begin(), names.end(), [&i]() { return std::to_string(i++); });
	}
	else {
		names = globalNames;
	}

	std::unordered_map<std::string, int> names2globalIds;

	for (size_t i = 0; i < names.size(); ++i) {
		names2globalIds[names[i]] = i;
	}

	std::vector<int> lowestClusterMembers(n_clusters, std::numeric_limits<int>::max());

	for (int local_id = 0; local_id < (int)assignments.size(); ++local_id) {
		// translate element ids
		int global_id = local2global[local_id];
		int cluster_id = assignments[local_id];

		if (global_id < lowestClusterMembers[cluster_id]) {
			lowestClusterMembers[cluster_id] = global_id;
		}
	}

	std::vector<std::string> representatives(names.size());

	int n_singletons = 0;
	for (size_t i = 0; i < names.size(); ++i) {

		int local_id = get_local_id(i);
		if (local_id == -1) {
			// not in matrix - own representative
			representatives[i] = names[i];
			++n_singletons;
			//	LOG_DEBUG << names[i] << endl;
		}
		else {
			int cluster_id = assignments[local_id];
			int lowest_member = lowestClusterMembers[cluster_id];
			representatives[i] = names[lowest_member];
		}
	}

	ofs << "object" << separator << "cluster" << std::endl;

	for (size_t i = 0; i < names.size(); ++i) {
		ofs << names[i] << separator << representatives[i] << std::endl;
	}

	return n_clusters + n_singletons;
}


/*********************************************************************************************************************/
template <class Distance>
int GraphNumbered<Distance>::saveAssignments(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator) const {

	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;
	int singleton_id = n_clusters;

	ofs << "object" << separator << "cluster" << std::endl;

	if (globalNames.empty()) {
		for (size_t local_id = 0; local_id < local2global.size(); ++local_id) {
			ofs << local2global[local_id] << "," << assignments[local_id] << std::endl;
		}
	}
	else {
		std::vector<int> globalAssignments(globalNames.size());


		for (size_t i = 0; i < globalNames.size(); ++i) {

			int local_id = get_local_id(i);
			if (local_id == -1) {
				// not in matrix 
				globalAssignments[i] = singleton_id++;
			}
			else {
				globalAssignments[i] = assignments[local_id];
			}
		}

		for (size_t i = 0; i < globalNames.size(); ++i) {
			ofs << globalNames[i] << separator << globalAssignments[i] << std::endl;
		}
	}

	return singleton_id; // return total number of clusters
}

