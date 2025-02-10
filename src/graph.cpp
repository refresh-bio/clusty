// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "graph.h"
#include "log.h"
#include "io.h"

#include <iostream>
#include <limits>
#include <unordered_map>
#include <algorithm>
#include <iterator>

using namespace std;

/*********************************************************************************************************************/
void Graph::sortClustersBySize(
	const std::vector<int>& assignments,
	std::vector<int>& old2new) const {

	// get number of clusters
	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;
	
	// calculate cluster sizes
	std::vector<std::pair<int, int>> clusters_n_sizes(n_clusters);
	int i = 0;
	std::generate(clusters_n_sizes.begin(), clusters_n_sizes.end(), [&i]() { return std::make_pair(i++, 0); });
	for (auto a : assignments) {
		++clusters_n_sizes[a].second;
	}

	// sort clusters decreasingly by size and establish mapping between cluster ids
	std::stable_sort(clusters_n_sizes.begin(), clusters_n_sizes.end(), [](const auto& p, const auto& q) { return p.second > q.second; });
	old2new.resize(n_clusters);
	for (int i = 0; i < n_clusters; ++i) {
		old2new[clusters_n_sizes[i].first] = i;
	}

}

/*********************************************************************************************************************/
void Graph::processHeader(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	const std::map<std::string, ColumnFilter>& columns2filters) {

	std::string line;
	std::getline(ifs, line);

	std::replace(line.begin(), line.end(), ',', ' ');
	std::istringstream iss(line);
	std::vector<std::string> columns;
	std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(columns));

	if (columns.size() < 3) {
		throw std::runtime_error("Error loading distances: at least three columns are required");
	}

	// default column identifiers
	int col_ids[]{ 0, 1 };
	int col_distance = 2; // by default use 3rd column as the one with distance 

	int n_columns = (int)columns.size();

	if (!idColumns.first.empty()) {
		col_ids[0] = std::find(columns.begin(), columns.end(), idColumns.first) - columns.begin();
		col_ids[1] = std::find(columns.begin(), columns.end(), idColumns.second) - columns.begin();

		if (col_ids[0] == n_columns || col_ids[1] == n_columns) {
			throw std::runtime_error("Error loading distances: id columns not found");
		}

		if (col_ids[0] > col_ids[1]) {
			std::swap(col_ids[0], col_ids[1]);
		}
	}

	if (!distanceColumn.empty()) {
		col_distance = std::find(columns.begin(), columns.end(), distanceColumn) - columns.begin();
		if (col_distance == n_columns) {
			throw std::runtime_error("Error loading distances: " + distanceColumn + " column not found");
		}
	}

	this->filters.resize(columns.size());

	for (const auto& f : columns2filters) {
		int col = std::find(columns.begin(), columns.end(), f.first) - columns.begin();
		if (col == n_columns) {
			throw std::runtime_error("Error loading distances: " + f.first + " column not found");
		}
		this->filters[col] = f.second;
		this->filters[col].enabled = true;
	}

	this->distanceColumnId = col_distance;
	this->sequenceColumnIds[0] = col_ids[0];
	this->sequenceColumnIds[1] = col_ids[1];
}

