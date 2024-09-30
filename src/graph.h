// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "distances.h"

#include <vector>
#include <limits>
#include <fstream>
#include <cstring>
#include <map>

// *******************************************************************************************/
struct ColumnFilter {
	double min{ std::numeric_limits<double>::lowest() };
	double max{ std::numeric_limits<double>::max() };
	bool enabled{ false };
};

// *******************************************************************************************/
class Graph {

public:
	virtual IMatrix& getMatrix() = 0;
	
	virtual size_t getNumVertices() const = 0;

	virtual size_t getNumInputVertices() const = 0;

	virtual size_t getNumEdges() const = 0;
		
	virtual size_t load(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		distance_transformation_t transform,
		const std::map<std::string, ColumnFilter>& columns2filters) = 0;

	virtual int saveAssignments(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator) const = 0;

	virtual int saveRepresentatives(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator) const = 0;

	virtual void reorderObjects(
		const std::vector<std::string>& externalNames,
		std::vector<int>& objects) const = 0;

	virtual void print(std::ostream& out) const = 0;

	virtual ~Graph() {}

protected:

	void processHeader(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		const std::map<std::string, ColumnFilter>& columns2filters,
		int idColumnsOut[2],
		int& distanceColumnOut,
		std::vector<ColumnFilter>& filters);
};