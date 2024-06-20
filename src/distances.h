#pragma once

#include <vector>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <map>
#include <iterator>
#include <functional>

#include "conversion.h"

using distance_transformation_t = std::function<double(double)>;

struct ColumnFilter {
	double min{ std::numeric_limits<double>::lowest() };
	double max{ std::numeric_limits<double>::max() };
	bool enabled{ false };
};

struct dist_t {
	static const double MAX;
	
	union {
		uint64_t ids;  ///< index (id)
		struct {
			uint32_t lo;   ///< numer wiersza
			uint32_t hi;   ///< numer kolumny
		} s;
	} u ;
	
	double d;  ///< distance   std::numeric_limits<double>::max() 

	dist_t() :  d(std::numeric_limits<double>::max()) { u.ids = 0; }
	dist_t(uint64_t ids, double d) :  d(d) { u.ids = ids; }
	dist_t(uint64_t i, uint64_t j, double d) :  d(d) { u.s.lo = i;  u.s.hi = j;}

	dist_t(const std::pair<uint64_t, double>& rhs) :  d(rhs.second) { u.ids = rhs.first;}

	static uint64_t pack(uint64_t i, uint64_t j) {
		if (i >= j) {
			std::swap(i, j);
		}
		return (i << 32ULL) | j;
	}

	static void unpack(uint64_t packed_ids, uint64_t& lo, uint64_t& hi) {
		hi = packed_ids & 0xffffffff;
		lo = packed_ids >> 32ULL;
	}

	// this is to preserve consistency with similarity variant
	// - increasingly by distance
	// - decreasingly by id
	bool operator<(const dist_t& rhs) const {
		return (d == rhs.d)
			? (u.ids > rhs.u.ids)
			: (d < rhs.d);
	}

	bool operator<=(const dist_t& rhs) const {
		return (d == rhs.d)
			? (u.ids >= rhs.u.ids)
			: (d <= rhs.d);
	}
};



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

class StringEqual {
	
public:
	StringEqual() {}

	bool operator()(const char* a, const char* b) const {
		return (std::strcmp(a, b) == 0);
	}
};


class SparseMatrix
{
protected:
	/** a vector of distances */
	std::vector<dist_t> distances;

	/** Each row holds a pointer to the first distance in 
	 *  the row. */
	std::vector<dist_t*> rows; /// wskaźniki na kolejne wiersze, ostatni pokazuje adres za kolekcją


public:

	SparseMatrix() {}
	
	virtual ~SparseMatrix() {}

	size_t num_objects() const { return rows.size() - 1; }

	size_t num_elements() const { return distances.size(); }
	void reserve(size_t count) {
		distances.reserve(count);
	}

	virtual size_t num_input_objects() const = 0;

	const dist_t* begin(int row_id) const { return rows[row_id]; }
	const dist_t* end(int row_id) const { return rows[row_id + 1]; }

   std::vector<std::tuple<std::size_t, std::size_t, double>> get_distances_of_objects() const 
   {
      std::vector<std::tuple<std::size_t, std::size_t, double>> objects;
      objects.reserve (distances.size());
      for (const auto & item : distances)
      {
         objects.emplace_back(item.u.s.lo, item.u.s.hi, item.d);
      }

      return objects;
   }
	const std::vector<dist_t> get_distances() const { return distances; }

	std::vector<dist_t> get_distances() { return distances; }

	size_t get_num_neighbours(int i) { return rows[i + 1] - rows[i]; }


	dist_t get(uint64_t i, uint64_t j) const 
	{
		dist_t v {i, j, dist_t::MAX };
		auto it = std::lower_bound(rows[i], rows[i + 1], v, [](const dist_t& a, const dist_t& b) { return a.u.s.hi < b.u.s.hi; });

		return (it == rows[i + 1] || it->u.s.hi != j ) ? v : *it;
	}

	void extract_row(uint32_t row, uint32_t count, dist_t* out) const 
	{
		dist_t* cur = rows[row];
		dist_t* end = rows[row + 1];
		
		for (uint32_t i = 0; i < count; ++i) 
		{
			if (i == row) {
				out[i].u.s.lo = out[i].u.s.hi = row;
				out[i].d = 0.0; // diagonal element
			}
			else if (cur == end || i < cur->u.s.hi) 
			{
				// before current or after end
				out[i].u.s.lo = row;
				out[i].u.s.hi = i;
				out[i].d = std::numeric_limits<double>::max();
			}
			else 
			{
				out[i] = *cur;
				++cur;
			}
		}
	}

	virtual size_t load(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		distance_transformation_t transform,
		std::map<std::string, ColumnFilter>& columns2filters) = 0;

	virtual int saveAssignments(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator) = 0;

	virtual int saveRepresentatives(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator) = 0;

	virtual void print(std::ostream& out) = 0;

protected:

	void processHeader(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		std::map<std::string, ColumnFilter>& columns2filters,
		int idColumnsOut[2],
		int& distanceColumnOut,
		std::vector<ColumnFilter>& filters);
	
};



class SparseMatrixNamed : public SparseMatrix {
	
	using ids_pair_t = std::pair<int, int>;
	
	std::unordered_map<const char*, ids_pair_t, StringHasher, StringEqual> names2ids;
	
	std::vector<const char*> ids2names;

	char* namesBuffer{ nullptr };

public:
	~SparseMatrixNamed() {
		delete[] namesBuffer;
	}

	size_t num_input_objects() const override { return names2ids.size(); }

	int get_id(const char* name) const {
		auto it = names2ids.find((char*)name);
		if (it == names2ids.end()) {
			return -1;
		}
		else {
			return it->second.first;
		}
	}

	const char* get_name(int id) const { return ids2names[id]; }

	size_t load(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		distance_transformation_t transform,
		std::map<std::string, ColumnFilter>& columns2filters) override;

	int saveAssignments(
		std::ofstream& ofs,
		const std::vector<std::string>& globalNames,
		const std::vector<int>& assignments,
		char separator) override;

	int saveRepresentatives(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator) override;

	void print(std::ostream& out) override;

};


class SparseMatrixNumbered : public SparseMatrix {
	
	std::vector<int> global2local;
	std::vector<int> local2global;

public:

	int get_local_id(int global_id) {
		if ((size_t)global_id > global2local.size() - 1) {
			return -1;
		}
		else {
			return global2local[global_id];
		}
	}
	
	size_t num_input_objects() const override { return global2local.size(); }

	size_t load(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		distance_transformation_t transform,
		std::map<std::string, ColumnFilter>& columns2filters) override;

	int saveAssignments(
		std::ofstream& ofs,
		const std::vector<std::string>& globalNames,
		const std::vector<int>& assignments,
		char separator);

	int saveRepresentatives(
		std::ofstream& ofs,
		const std::vector<std::string>& externalNames,
		const std::vector<int>& assignments,
		char separator);

	void print(std::ostream& out) override {}
};
