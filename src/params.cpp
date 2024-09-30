// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "params.h"
#include "distances.h"
#include "log.h"
#include "leiden.h"

#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

void Params::printUsage() const {
	LOG_NORMAL << "Usage:" << endl
		<< "clusty [options] <distances> <assignments>" << endl << endl
		<< "Parameters:" << endl
		<< "  <distances> - input TSV/CSV table with pairwise distances" << endl
		<< "  <assignments> - output TSV/CSV table with assignments" << endl << endl
		<< "Options:" << endl
		<< "  " + PARAM_FILE_OBJECTS + " <string> - optional TSV/CSV file with object names in the first column sorted decreasingly w.r.t. representativness" << endl
		<< "  " + PARAM_ALGO + " <single | complete | uclust | set-cover | cd-hit | leiden> - clustering algorithm:" << endl
		<< "    * single     - single linkage (connected component, MMseqs mode 1)" << endl
		<< "    * complete   - complete linkage" << endl
		<< "    * uclust     - UCLUST" << endl
		<< "    * set-cover  - greedy set cover (MMseqs mode 0)" << endl
		<< "    * cd-hit     - CD-HIT (greedy incremental; MMseqs mode 2)" << endl
		<< "    * leiden     - Leiden algorithm" << endl << endl

		<< "  " + PARAM_ID_COLUMNS + " <column-name1> <column-name2> - names of columns with sequence identifiers (default: two first columns)" << endl
		<< "  " + PARAM_DISTANCE_COLUMN + " <column-name> - name of the column with pairwise distances (or similarities; default: third column)" << endl
		<< "  " + FLAG_SIMILARITY + " - use similarity instead of distances (has to be in [0,1] interval; default: false)" << endl
		<< "  " + FLAG_PERCENT_SIMILARITY + " - use percent similarity (has to be in [0,100] interval; default: false)" << endl
		<< "  " + PARAM_MIN + " <column-name> <real-threshold> - accept pairwise connections with values greater or equal given threshold in a specified column" << endl
		<< "  " + PARAM_MAX + " <column-name> <real-threshold> - accept pairwise connections with values lower or equal given threshold in a specified column" << endl
		<< "  " + FLAG_NUMERIC_IDS + " - use when sequences in the distances file are represented by numbers (can be mapped to string ids by the object file)" << endl
		<< "  " + FLAG_OUT_REPRESENTATIVES + " - output a representative object for each cluster instead of a cluster numerical identifier (default: " << std::boolalpha << outputRepresentatives << ")" << endl
		<< "  " + FLAG_OUT_CSV + " - output a CSV table instead of a default TSV (default: " << std::boolalpha << outputCSV << ")"

#ifndef NO_LEIDEN
		<< endl << endl
		<< "  " + PARAM_LEIDEN_RESOLUTION + " - resolution parameter for Leiden algorithm (default: " << leidenParams.resolution << ")" << endl
		<< "  " + PARAM_LEIDEN_BETA + " - beta parameter for Leiden algorithm (default: " << leidenParams.beta << ")" << endl
		<< "  " + PARAM_LEIDEN_RESOLUTION + " - number of interations for Leiden algorithm (default: " << leidenParams.numIterations << ")"
#endif
		<< endl << endl;
}


bool Params::parse(int argc, char** argv) {

	vector<string> args;
	for (int i = 1; i < argc; ++i) {
		args.emplace_back(argv[i]);
	}

	if (args.size() < 2) {
		return false;
	}

	std::string tmp;
	findOption(args, PARAM_ALGO, tmp);
	if (tmp.length()) {
		algo = str2algo(tmp);
	}

	findOption(args, PARAM_FILE_OBJECTS, objectsFile);

	findOption(args, PARAM_ID_COLUMNS, idColumns.first, idColumns.second);
	numericIds = findSwitch(args, FLAG_NUMERIC_IDS);

	findOption(args, PARAM_DISTANCE_COLUMN, distanceColumn);

	bool use_similarity = findSwitch(args, FLAG_SIMILARITY);
	bool use_percent_similarity = findSwitch(args, FLAG_PERCENT_SIMILARITY);

	if (use_percent_similarity) {
		distanceSpecification = DistanceSpecification::PercentSimilarity;
	}
	else if (use_similarity) {
		distanceSpecification = DistanceSpecification::Similarity;
	}

	string column;
	double value;
	while (findOption(args, PARAM_MIN, column, value)) {
		columns2filters[column].min = std::max(value, columns2filters[column].min);
	}
	while (findOption(args, PARAM_MAX, column, value)) {
		columns2filters[column].max = std::min(value, columns2filters[column].max);
	}

	outputRepresentatives = findSwitch(args, FLAG_OUT_REPRESENTATIVES);
	outputCSV = findSwitch(args, FLAG_OUT_CSV);

	// leiden parameters
	findOption(args, PARAM_LEIDEN_RESOLUTION, leidenParams.resolution);
	findOption(args, PARAM_LEIDEN_BETA, leidenParams.beta);
	findOption(args, PARAM_LEIDEN_ITERATIONS, leidenParams.numIterations);

	verbose = findSwitch(args, FLAG_VERBOSE);

	if (args.size() == 2) {
		distancesFile = args[0];
		output = args[1];
		return true;
	}

	return false;
}

