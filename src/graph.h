// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "distances.h"

#include <vector>
#include <limits>
#include <fstream>
#include <cstring>
#include <map>
#include <tuple>

// *******************************************************************************************/
struct ColumnFilter {
	double min{ std::numeric_limits<double>::lowest() };
	double max{ std::numeric_limits<double>::max() };
	bool enabled{ false };
};

// *******************************************************************************************/
class IEdgesCollection {
public:
	virtual ~IEdgesCollection() {}

	virtual void clear() = 0;
};

// *******************************************************************************************/
template <class edge_label_t>
class EdgesCollection : public IEdgesCollection {
public:
	using edge_t = std::pair<edge_label_t[2], double>;

	EdgesCollection(size_t preallocSize) {
		data.reserve(preallocSize); 
	}

	void clear() override { data.clear(); }

	std::vector<edge_t> data;
	
	int maxEdge{ 0 };
};


// *******************************************************************************************/
class Graph {

protected:
	int numThreads; 

	int sequenceColumnIds[2]{ 0, 1 };

	int distanceColumnId{ 2 };

	std::vector<ColumnFilter> filters;

public:
	static bool isSeparator(char c) { return c == ',' || c == '\t' || c == '\r' || c == '\n'; }
	static bool isNewline(char c) { return c == '\r' || c == '\n'; }

	Graph(int numThreads) : numThreads(std::max(4, numThreads)) {}; // at least three - loader,parser,updater
		 
	virtual ~Graph() {}

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
		const std::vector<std::string_view>& externalNames,
		const std::vector<int>& assignments,
		char separator,
		bool useRepresentatives) const = 0;

	virtual void reorderObjects(
		const std::vector<std::string_view>& externalNames,
		std::vector<int>& objects) const = 0;

	virtual void print(std::ostream& out) const = 0;



protected:

	void processHeader(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		const std::map<std::string, ColumnFilter>& columns2filters);


	void sortClustersBySize(
		const std::vector<int>& assignments,
		std::vector<int>& old2new
	) const;


	template <class object_t, class... Ts>
	void fillRepresentatives(
		const std::vector<std::tuple<object_t, int, Ts...>>& objects_and_clusters,
		std::vector<std::tuple<object_t, object_t>>& objects_and_representatives) const;
	
	
};



template <class object_t, class... Ts>
void Graph::fillRepresentatives(
	const std::vector<std::tuple<object_t, int, Ts...>> & objects_and_clusters,
	std::vector<std::tuple<object_t, object_t>>& objects_and_representatives) const {

	auto representative = std::get<0>(objects_and_clusters.front());
	int cluster = std::get<1>(objects_and_clusters.front());

	objects_and_representatives.resize(objects_and_clusters.size());

	for (int i = 0; i < (int)objects_and_clusters.size(); ++i) {
		const auto& in = objects_and_clusters[i];
		auto& out = objects_and_representatives[i];
	
		// update representative
		if (std::get<1>(in) != cluster) {
			cluster = std::get<1>(in);
			representative = std::get<0>(in);
		}
		
		std::get<0>(out) = std::get<0>(in);
		std::get<1>(out) = representative;
	}
}