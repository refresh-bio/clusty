// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "graph.h"
#include "leiden.h"

#include <string>
#include <stdexcept>
#include <sstream>
#include <map>

enum class Algo {
	SingleLinkage,
	CompleteLinkage,
	UClust,
	SetCover,
	CdHit,
	Leiden
};

enum class DistanceSpecification {
	Distance,
	Similarity,
	PercentSimilarity
};


class Params {
	const std::string PARAM_ALGO{ "--algo" };

	const std::string PARAM_FILE_OBJECTS{ "--objects-file" };

	const std::string PARAM_ID_COLUMNS{ "--id-cols" };
	const std::string FLAG_NUMERIC_IDS{ "--numeric-ids" };
	const std::string PARAM_DISTANCE_COLUMN{ "--distance-col" };
	const std::string FLAG_SIMILARITY{ "--similarity" };
	const std::string FLAG_PERCENT_SIMILARITY{ "--percent-similarity" };

	const std::string PARAM_MAX{ "--max" };
	const std::string PARAM_MIN{ "--min" };

	const std::string FLAG_OUT_REPRESENTATIVES{ "--out-representatives" };
	const std::string FLAG_OUT_CSV{ "--out-csv" };

	const std::string PARAM_LEIDEN_RESOLUTION{ "--leiden-resolution" };
	const std::string PARAM_LEIDEN_BETA{ "--leiden-beta" };
	const std::string PARAM_LEIDEN_ITERATIONS{ "--leiden-iterations" };

	const std::string FLAG_VERBOSE{ "-v" };
	const std::string FLAG_VERSION{ "--version" };

	const std::string PARAM_NUM_THREADS{ "-t" };

public:
	enum Status {
		Correct,
		Incorrect,
		ShowVersion
	};
	
	static Algo str2algo(const std::string& str)
	{
		if (str == "single") { return Algo::SingleLinkage; }
		else if (str == "complete") { return Algo::CompleteLinkage; }
		else if (str == "uclust") { return Algo::UClust; }
		else if (str == "set-cover") { return Algo::SetCover; }
		else if (str == "cd-hit") { return Algo::CdHit; }
		else if (str == "leiden") { return Algo::Leiden; }

		else { throw std::runtime_error("Unkown clustering algorithm"); }
	}

	static std::string algo2str(Algo algo) {
		switch (algo) {
		case Algo::SingleLinkage: return "single";
		case Algo::CompleteLinkage: return "complete";
		case Algo::UClust: return "uclust";
		case Algo::SetCover: return "set-cover";
		case Algo::Leiden: return "leiden";
		case Algo::CdHit: return "cd-hit";
		default: throw std::runtime_error("Unkown clustering algorithm");
		}
	}


public:

	std::string objectsFile;
	std::string distancesFile;
	std::string output;

	Algo algo{ Algo::SingleLinkage };

	std::pair<std::string, std::string> idColumns;
	bool numericIds{ false };

	std::string distanceColumn;
	DistanceSpecification distanceSpecification{ DistanceSpecification::Distance };

	std::map<std::string, ColumnFilter> columns2filters;
	bool outputRepresentatives{ false };
	bool outputCSV{ false };

	LeidenParams leidenParams;

	bool verbose{ false };

	int numThreads{ 4 };

	void printUsage() const;
	Status parse(int argc, char** argv);

	bool findSwitch(std::vector<std::string>& params, const std::string& name) {
		auto it = find(params.begin(), params.end(), name); // verbose mode
		if (it != params.end()) {
			params.erase(it);
			return true;
		}

		return false;
	}

	template <typename T>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& v) {
		if (params.size() < 2) {
			return false;
		}

		auto stop = std::prev(params.end());
		auto it = find(params.begin(), stop, name); // verbose mode
		if (it != stop) {
			std::istringstream iss(*std::next(it));
			if (iss >> v) {
				params.erase(it, it + 2);
				return true;
			}
		}
		return false;
	}

	template <typename T, typename U>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& value1, U& value2) {
		if (params.size() < 3) {
			return false;
		}

		auto stop = std::prev(params.end(), 2);
		auto it = find(params.begin(), stop, name); // verbose mode
		if (it != stop) {
			if (std::istringstream(*std::next(it)) >> value1
				&& std::istringstream(*std::next(it, 2)) >> value2) {
				params.erase(it, it + 3);
				return true;
			}
		}
		return false;
	}
};
