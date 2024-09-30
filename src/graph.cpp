// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "graph.h"
#include "log.h"

#include <iostream>
#include <limits>
#include <unordered_map>
#include <algorithm>
#include <iterator>

using namespace std;

/*********************************************************************************************************************/
void Graph::processHeader(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	const std::map<std::string, ColumnFilter>& columns2filters,
	int idColumnsOut[2],
	int& distanceColumnOut,
	std::vector<ColumnFilter>& filters) {

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

	filters.resize(columns.size());

	for (const auto& f : columns2filters) {
		int col = std::find(columns.begin(), columns.end(), f.first) - columns.begin();
		if (col == n_columns) {
			throw std::runtime_error("Error loading distances: " + f.first + " column not found");
		}
		filters[col] = f.second;
		filters[col].enabled = true;
	}

	idColumnsOut[0] = col_ids[0];
	idColumnsOut[1] = col_ids[1];
	distanceColumnOut = col_distance;
}

