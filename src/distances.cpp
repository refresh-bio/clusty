#include "distances.h"
#include "log.h"

#include <iostream>
#include <limits>
#include <unordered_map>

using namespace std;

const double dist_t::MAX = std::numeric_limits<double>::max();

/*********************************************************************************************************************/
void SparseMatrix::processHeader(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	std::map<std::string, ColumnFilter>& columns2filters,
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

/*********************************************************************************************************************/
size_t SparseMatrixNamed::load(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	distance_transformation_t transform,
	std::map<std::string, ColumnFilter>& columns2filters) {

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

	std::vector<int> counts;
	std::vector<dist_t> tmp_dists;
	tmp_dists.reserve(128LL << 20); // for 128M distances

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
			decltype(names2ids)::iterator its[2];
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
					counts.push_back(0);
					it->second.first = ids2names.size() - 1;
				}
			}

			uint32_t i = std::min(its[0]->second.first, its[1]->second.first);
			uint32_t j = std::max(its[0]->second.first, its[1]->second.first);

			// omit diagonal elements - they are assumed to have 0 distance
			if (i == j) {
				continue;
			}

			++counts[i];
			++counts[j];

			tmp_dists.emplace_back(i, j, d);
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


	distances.resize(tmp_dists.size() * 2);
	rows.resize(counts.size());
	int cumulated = 0;
	for (size_t i = 0; i < rows.size(); ++i) {
		rows[i] = distances.data() + cumulated;
		cumulated += counts[i];
	}

	struct row_info {
		int n_filled{ 0 };
		int last_id{ -1 };
	};

	std::vector <row_info> rows_info(rows.size());

	// second pass - put distances in the final structure
	for (const dist_t& dist : tmp_dists) {

		uint32_t i = dist.u.s.lo;
		uint32_t j = dist.u.s.hi;
		double d = dist.d;

		rows[i][rows_info[i].n_filled] = dist;
		++rows_info[i].n_filled;
		rows_info[i].last_id = j;

		rows[j][rows_info[j].n_filled] = dist_t{ j,i,d };
		++rows_info[j].n_filled;
		rows_info[j].last_id = i;
	}

	auto end = distances.data() + distances.size();
	rows.push_back(end);

	// if neccessary, sort distances in rows according to the second id
	dist_t* curBegin = rows[0];

	for (size_t i = 0; i < rows.size() - 1; ++i) {
		std::sort(curBegin, rows[i + 1], [](const dist_t& a, const dist_t& b) { return a.u.ids < b.u.ids; });
		auto newEnd = std::unique(curBegin, rows[i + 1], [](const dist_t& a, const dist_t& b) { return a.u.ids == b.u.ids; });
		
		if (rows[i] != curBegin) {
			newEnd = std::copy(curBegin, newEnd, rows[i]);
		}
		
		curBegin = rows[i + 1];
		rows[i + 1] = newEnd;
	}

	size_t newSize = rows.back() - rows.front();
	distances.erase(distances.begin() + newSize, distances.end());

	delete[] buf;

	// debug stuff
	//std::ofstream dbg("debug.log");
	//print(dbg);

	return n_total_distances;
}

/*********************************************************************************************************************/
int SparseMatrixNamed::saveAssignments(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator) {

	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;

	vector<string> names;
	if (globalNames.empty()) {
		names.resize(names2ids.size());
		vector<pair<const char*, pair<int,int>>> entries(names2ids.size());
		copy(names2ids.begin(), names2ids.end(), entries.begin());
		sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.second.second < b.second.second; });

		transform(entries.begin(), entries.end(), names.begin(), [](const auto& entry) { return string(entry.first); });
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

	ofs << "object" << separator << "cluster" << endl;

	for (size_t i = 0; i < names.size(); ++i) {
		ofs << names[i] << separator << globalAssignments[i] << endl;
	}

	return singleton_id;
}

/*********************************************************************************************************************/
int SparseMatrixNamed::saveRepresentatives(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator
) {
	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;
	
	std::vector<std::string> names;
	if (globalNames.empty()) {
		names.resize(names2ids.size());
		vector<pair<const char*, pair<int, int>>> entries(names2ids.size());
		copy(names2ids.begin(), names2ids.end(), entries.begin());
		sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.second.second < b.second.second; });

		transform(entries.begin(), entries.end(), names.begin(), [](const auto& entry) { return string(entry.first); });
	}
	else {
		names = globalNames;
	}
	
	std::unordered_map<string, int> names2globalIds;

	for (size_t i = 0; i < names.size(); ++i) {
		names2globalIds[names[i]] = i;
	}

	std::vector<int> lowestClusterMembers(n_clusters, std::numeric_limits<int>::max());

	for (size_t i = 0; i < assignments.size(); ++i) {
		// translate element ids
		string name = get_name(i);
		int global_id = names2globalIds[name];
		int cluster_id = assignments[i];

		if (global_id < lowestClusterMembers[cluster_id]) {
			lowestClusterMembers[cluster_id] = global_id;
		}
	}

	std::vector<string*> representatives(names.size());

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

	ofs << "object" << separator << "cluster" << endl;

	for (size_t i = 0; i < names.size(); ++i) {
		ofs << names[i] << separator << *representatives[i] << endl;
	}

	return n_clusters + n_singletons;
}


/*********************************************************************************************************************/
size_t SparseMatrixNumbered::load(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	distance_transformation_t transform,
	std::map<std::string, ColumnFilter>& columns2filters) {

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

	std::vector<int> counts;

	counts.reserve(8LL << 20); // assume space for 8M objects
	global2local.reserve(8LL << 20);
	local2global.reserve(8LL << 20);

	std::vector<dist_t> tmp_dists;
	tmp_dists.reserve(128LL << 20); // for 128M distances

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

			// resize counts vector
			if (j + 1 > counts.size()) {
				counts.resize(j + 1);
			}

			++counts[i];
			++counts[j];

			tmp_dists.emplace_back(i, j, d);
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


	distances.resize(tmp_dists.size() * 2);
	rows.resize(counts.size());
	int cumulated = 0;
	for (size_t i = 0; i < rows.size(); ++i) {
		rows[i] = distances.data() + cumulated;
		cumulated += counts[i];
	}

	struct row_info {
		int n_filled{ 0 };
		int last_id{ -1 };
	};

	std::vector <row_info> rows_info(rows.size());

	// second pass - put distances in the final structure
	for (const dist_t& dist : tmp_dists) {

		uint32_t i = dist.u.s.lo;
		uint32_t j = dist.u.s.hi;
		double d = dist.d;

		rows[i][rows_info[i].n_filled] = dist;
		++rows_info[i].n_filled;
		rows_info[i].last_id = j;

		rows[j][rows_info[j].n_filled] = dist_t{ j,i,d };
		++rows_info[j].n_filled;
		rows_info[j].last_id = i;
	}

	auto end = distances.data() + distances.size();
	rows.push_back(end);

	// if neccessary, sort distances in rows according to the second id
	dist_t* curBegin = rows[0];

	for (size_t i = 0; i < rows.size() - 1; ++i) {
		std::sort(curBegin, rows[i + 1], [](const dist_t& a, const dist_t& b) { return (a.u.ids == b.u.ids) ? (a.d < b.d) : (a.u.ids < b.u.ids); });
		auto newEnd = std::unique(curBegin, rows[i + 1], [](const dist_t& a, const dist_t& b) { return a.u.ids == b.u.ids; });

		if (rows[i] != curBegin) {
			newEnd = std::copy(curBegin, newEnd, rows[i]);
		}

		curBegin = rows[i + 1];
		rows[i + 1] = newEnd;
	}

	size_t newSize = rows.back() - rows.front();
	distances.erase(distances.begin() + newSize, distances.end());

	delete[] buf;

	return n_total_distances;
}


/*********************************************************************************************************************/
int SparseMatrixNumbered::saveAssignments(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator) {

	int n_clusters = *std::max_element(assignments.begin(), assignments.end()) + 1;
	int singleton_id = n_clusters;

	ofs << "object" << separator << "cluster" << endl;

	if (globalNames.empty()) {
		for (size_t local_id = 0; local_id < local2global.size(); ++local_id) {
			ofs << local2global[local_id] << "," << assignments[local_id] << endl;
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
			ofs << globalNames[i] << separator << globalAssignments[i] << endl;
		}
	}

	return singleton_id; // return total number of clusters
}


/*********************************************************************************************************************/
int SparseMatrixNumbered::saveRepresentatives(
	std::ofstream& ofs,
	const std::vector<std::string>& globalNames,
	const std::vector<int>& assignments,
	char separator
) {

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

	std::unordered_map<string, int> names2globalIds;

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

	std::vector<string> representatives(names.size());

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

	ofs << "object" << separator << "cluster" << endl;

	for (size_t i = 0; i < names.size(); ++i) {
		ofs << names[i] << separator << representatives[i] << endl;
	}

	return n_clusters + n_singletons;
}

/*********************************************************************************************************************/
void SparseMatrixNamed::print(std::ostream& out) {

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

		int i = names2ids[name].first;
		std::vector<dist_t> row(rows[i], rows[i + 1]);

		std::sort(row.begin(), row.end(), [this](const dist_t& a, const dist_t& b) {
			return strcmp(ids2names[a.u.s.hi], ids2names[b.u.s.hi]) < 0;
			});

		for (auto& p : row) {
			out << ids2names[p.u.s.lo] << "," << ids2names[p.u.s.hi] << "," << std::fixed << p.d << std::endl;
		}
	}
}

