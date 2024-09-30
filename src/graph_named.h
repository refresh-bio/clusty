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

#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>



// *******************************************************************************************
class StringHasher {
	std::hash<char> hasher;
public:
	StringHasher() {}

	size_t operator()(const char* s) const {
		size_t hs = hasher(*s);
		++s;
		while (*s) {
			hs ^= hasher(*s);
			++s;
		}
		return hs;
	}
};

// *******************************************************************************************
class StringEqual {

public:
	StringEqual() {}

	bool operator()(const char* a, const char* b) const {
		return (std::strcmp(a, b) == 0);
	}
};

/*********************************************************************************************************************/
template <class Distance>
class GraphNamed : public Graph {

	using ids_pair_t = std::pair<int, int>;

	SparseMatrix<Distance> matrix;

	std::unordered_map<const char*, ids_pair_t, StringHasher, StringEqual> names2ids;

	std::vector<const char*> ids2names;

	char* namesBuffer{ nullptr };

public:
	~GraphNamed() {
		delete[] namesBuffer;
	}

	IMatrix& getMatrix() override { return matrix; }

	size_t getNumVertices() const override { return matrix.num_objects(); }

	size_t getNumInputVertices() const override { return this->names2ids.size(); }

	size_t getNumEdges() const override { return matrix.num_elements(); }

	void reorderObjects(
		const std::vector<std::string>& externalNames,
		std::vector<int>& objects) const override {

		int obj_id = 0;
		for (const std::string& name : externalNames) {
			int local_id = get_id(name.c_str());
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

	void print(std::ostream& out) const override;

protected:
	
	int get_id(const char* name) const {
		auto it = names2ids.find((char*)name);
		if (it == names2ids.end()) {
			return -1;
		}
		else {
			return it->second.first;
		}
	}
};



/*********************************************************************************************************************/
template <class Distance>
size_t GraphNamed<Distance>::load(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	distance_transformation_t transform,
	const std::map<std::string, ColumnFilter>& columns2filters) {

	size_t n_total_distances = 0;
	const size_t N = 512ULL << 20; // 128MB buffer
	char* buf = new char[N];
	char* buf_end = buf + N;

	// get header
	int col_ids[]{ 0, 1 };
	int col_distance = 2; // by default use 3rd column as the one with distance 
	std::vector<ColumnFilter> filters;
	processHeader(ifs, idColumns, distanceColumn, columns2filters, col_ids, col_distance, filters);
	int n_columns = filters.size();

	auto is_sep = [](char c) {return c == ',' || c == '\t' || c == '\r' || c == '\t'; };
	auto is_newline = [](char c) {return c == '\r' || c == '\n'; };

	namesBuffer = new char[1LL << 30]; // 1 GB buffer for names
	char* raw_ptr = namesBuffer;

	// assume space for 8M objects
	matrix.distances.reserve(8LL << 20);

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

		//LOG_DEBUG << "portion: " << n_read << ", offset: " << offset << ", carryOn: " << carryOn << endl;

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
			decltype(names2ids.begin()) its[2];
			double d = std::numeric_limits<double>::max();

			int k = 0;
			bool carryOn = true;

			for (int c = 0; c < n_columns; ++c) {
				char* q = std::find_if(p, line_end, is_sep); // support both tsv and csv files
				*q = 0;

				if (k < 2 && c == col_ids[k]) {
					char* name = p;

					// store name in hashtable
					its[k] = names2ids.find(name);

					if (its[k] == names2ids.end()) {

						char* localName = raw_ptr;
						while (*raw_ptr++ = *name++) {}
						auto it_and_flag = names2ids.insert({ localName, {-1, names2ids.size()} }); // -1 indicate singleton
						its[k] = it_and_flag.first;
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

			// translate seq names to numerical ids
			for (int k = 0; k < 2; ++k) {
				auto it = its[k];

				// if name not mapped to numerical ids
				if (it->second.first == -1) {
					ids2names.push_back(it->first);
					it->second.first = ids2names.size() - 1;
				}
			}

			uint32_t i = std::min(its[0]->second.first, its[1]->second.first);
			uint32_t j = std::max(its[0]->second.first, its[1]->second.first);

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

	return n_total_distances;
}



/*********************************************************************************************************************/
template <class Distance>
int GraphNamed<Distance>::saveRepresentatives(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator) const
{
	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;

	std::vector<std::string> names;
	if (globalNames.empty()) {
		names.resize(names2ids.size());
		std::vector<std::pair<const char*, std::pair<int, int>>> entries(names2ids.size());
		std::copy(names2ids.begin(), names2ids.end(), entries.begin());
		std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.second.second < b.second.second; });

		std::transform(entries.begin(), entries.end(), names.begin(), [](const auto& entry) { return std::string(entry.first); });
	}
	else {
		names = globalNames;
	}

	std::unordered_map<std::string, int> names2globalIds;

	for (size_t i = 0; i < names.size(); ++i) {
		names2globalIds[names[i]] = i;
	}

	std::vector<int> lowestClusterMembers(n_clusters, std::numeric_limits<int>::max());

	for (size_t i = 0; i < assignments.size(); ++i) {
		// translate element ids
		std::string name = ids2names[i];
		int global_id = names2globalIds[name];
		int cluster_id = assignments[i];

		if (global_id < lowestClusterMembers[cluster_id]) {
			lowestClusterMembers[cluster_id] = global_id;
		}
	}

	std::vector<std::string*> representatives(names.size());

	int n_singletons = 0;
	for (size_t i = 0; i < names.size(); ++i) {

		int id = get_id(names[i].c_str());
		if (id == -1) {
			// not in matrix - own representative
			representatives[i] = &names[i];
			++n_singletons;
			//	LOG_DEBUG << names[i] << endl;
		}
		else {
			int cluster_id = assignments[id];
			int lowest_member = lowestClusterMembers[cluster_id];
			representatives[i] = &names[lowest_member];
		}
	}

	ofs << "object" << separator << "cluster" << std::endl;

	for (size_t i = 0; i < names.size(); ++i) {
		ofs << names[i] << separator << *representatives[i] << std::endl;
	}

	return n_clusters + n_singletons;
}



/*********************************************************************************************************************/
template <class Distance>
int GraphNamed<Distance>::saveAssignments(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator) const {

	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;

	std::vector<std::string> names;
	if (globalNames.empty()) {
		names.resize(names2ids.size());
		std::vector<std::pair<const char*, std::pair<int, int>>> entries(names2ids.size());
		std::copy(names2ids.begin(), names2ids.end(), entries.begin());
		std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.second.second < b.second.second; });

		std::transform(entries.begin(), entries.end(), names.begin(), [](const auto& entry) { return std::string(entry.first); });
	}
	else {
		names = globalNames;
	}

	std::vector<int> globalAssignments(names.size());

	int singleton_id = n_clusters;

	for (size_t i = 0; i < names.size(); ++i) {

		int id = get_id(names[i].c_str());
		if (id == -1) {
			// not in matrix 
			globalAssignments[i] = singleton_id++;
		}
		else {
			globalAssignments[i] = assignments[id];
		}
	}

	ofs << "object" << separator << "cluster" << std::endl;

	for (size_t i = 0; i < names.size(); ++i) {
		ofs << names[i] << separator << globalAssignments[i] << std::endl;
	}

	return singleton_id;
}

/*********************************************************************************************************************/
template <class Distance>
void GraphNamed<Distance>::print(std::ostream& out) const {

	out.precision(10);

	std::vector<const char*> names(names2ids.size());
	int i = 0;
	for (auto q : names2ids) {
		names[i] = q.first;
		++i;
	}

	std::sort(names.begin(), names.end(), [](const char* a, const char* b) {
		return strcmp(a, b) < 0;
		});

	for (auto name : names) {

		int i = names2ids.at(name).first;
		const std::vector<Distance>& row = matrix.distances[i];

		for (auto& p : row) {
			out << ids2names[i] << "," << ids2names[p.get_id()] << "," << std::fixed << p.get_d() << std::endl;
		}
	}
}

