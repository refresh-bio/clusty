// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "console.h"
#include "log.h"

#include "graph_named.h"
#include "graph_numbered.h"
#include "sparse_matrix.h"

#define VAL(str) #str
#define TOSTRING(str) VAL(str)
#include "version.h"

#include <vector>
#include <memory>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

// *******************************************************************************************
bool Console::init(int argc, char** argv, Params& params) {
	Log::getInstance(Log::LEVEL_NORMAL).enable();

	LOG_NORMAL << "Clusty" << endl
		<< "  version " << VERSION
#ifdef GIT_COMMIT
		<< "-" << TOSTRING(GIT_COMMIT)
#endif
		<< " (" << DATE << ")" << endl << endl;

	if (!params.parse(argc, argv)) {
		params.printUsage();
		return false;
	}

	if (params.verbose) {
		Log::getInstance(Log::LEVEL_VERBOSE).enable();
	}

	return true;
}

// *******************************************************************************************
std::unique_ptr<Graph> Console::loadGraph(const Params& params) {
	
	unique_ptr<Graph> graph;

	if (params.numericIds) {
		if (needDistances(params.algo)) {
			graph = make_unique<GraphNumbered<dist_t>>();
		}
		else {
			graph = make_unique<GraphNumbered<mini_dist_t>>();
		}
	}
	else {
		if (needDistances(params.algo)) {
			graph = make_unique<GraphNamed<dist_t>>();
		}
		else {
			graph = make_unique<GraphNamed<mini_dist_t>>();
		}
	}

	LOG_NORMAL << "Loading pairwise distances from " << params.distancesFile << "... ";
	auto t = std::chrono::high_resolution_clock::now();

	vector<char> filebuf(128ULL << 20);  // 128MB buffer
	ifstream ifs;
	ifs.rdbuf()->pubsetbuf(filebuf.data(), filebuf.size());
	ifs.open(params.distancesFile, ios_base::binary);

	if (!ifs) {
		throw std::runtime_error("Unable to open distance file");
	}

	map<DistanceSpecification, distance_transformation_t> transforms{
		{ DistanceSpecification::Distance, [](double d) { return d; } },
		{ DistanceSpecification::Similarity,		[](double d) { return 1.0 - d; } },
		{ DistanceSpecification::PercentSimilarity, [](double d) { return 1.0 - d * 0.01; } },
	};

	size_t n_total_dists = graph->load(ifs, params.idColumns, params.distanceColumn,
		transforms[params.distanceSpecification], params.columns2filters);
	
	auto dt = std::chrono::high_resolution_clock::now() - t;

	ifs.close();
	LOG_NORMAL << endl
		<< "  input graph: " << graph->getNumInputVertices() << " nodes, " << n_total_dists << " edges" << endl
		<< "  filtered graph: " << graph->getNumVertices() << " nodes, " << graph->getNumEdges() << " edges" << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;

	return graph;
}

// *******************************************************************************************
void Console::loadObjects(
	const Params& params,
	const Graph& graph,
	std::vector<int>& objects,
	std::vector<std::string>& names) {

	names.clear();
	objects.resize(graph.getNumVertices());
	std::iota(objects.begin(), objects.end(), 0);
	ifstream ifs;

	if (!params.objectsFile.empty()) {
		LOG_NORMAL << "Loading objects from " << params.objectsFile << "... ";
		ifs.open(params.objectsFile);

		if (ifs) {
			auto t = std::chrono::high_resolution_clock::now();
			string token;
			token.reserve(16 << 10); // 16k buffer

			// omit header
			getline(ifs, token);
			auto is_sep = [](char c) {return c == ',' || c == '\t' || c == '\r' || c == '\t'; };

			while (getline(ifs, token)) {
				auto it = find_if(token.begin(), token.end(), is_sep);
				if (it != token.begin()) {
					names.emplace_back(token.begin(), it);
				}
			}

			ifs.close();
			LOG_NORMAL << endl;

			graph.reorderObjects(names, objects);

			auto dt = std::chrono::high_resolution_clock::now() - t;
			LOG_NORMAL << "  total objects: " << names.size() << endl
				<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
		}
		else {
			throw std::runtime_error("Unable to open objects file");
		}
	}
}


// *******************************************************************************************
void Console::doClustering(
	const Params& params,
	Graph& graph,
	const std::vector<int>& objects,
	std::vector<int>& assignments)
{
	assignments.clear();
	
	LOG_NORMAL << "Clustering (algorithm: " << Params::algo2str(params.algo) << ")... ";

	auto t = std::chrono::high_resolution_clock::now();
	double threshold = std::nexttoward(std::numeric_limits<double>::max(), 0.0);
	int n_clusters = 0;

	if (needDistances(params.algo)) {
		auto clustering = createClusteringAlgo<dist_t>(params);
		IMatrix& mat = graph.getMatrix();
		SparseMatrix<dist_t>& distances = static_cast<SparseMatrix<dist_t>&>(mat);
		n_clusters = (*clustering)(distances, objects, threshold, assignments);
	}
	else {
		auto clustering = createClusteringAlgo<mini_dist_t>(params);
		IMatrix& mat = graph.getMatrix();
		SparseMatrix<mini_dist_t>& distances = static_cast<SparseMatrix<mini_dist_t>&>(mat);
		n_clusters = (*clustering)(distances, objects, threshold, assignments);
	}

	auto dt = std::chrono::high_resolution_clock::now() - t;
	
	LOG_NORMAL << endl
		<< "  objects: " << graph.getNumVertices() << ", clusters: " << n_clusters << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
}

// *******************************************************************************************
void Console::saveAssignments(
	const Params& params,
	const Graph& graph,
	const std::vector<std::string>& names,
	const std::vector<int>& assignments) {

	LOG_NORMAL << "Saving clusters... ";
	auto t = std::chrono::high_resolution_clock::now();

	char sep = params.outputCSV ? ',' : '\t';

	ofstream ofs(params.output);
	int n_total_clusters = 0;
	if (params.outputRepresentatives) {
		n_total_clusters = graph.saveRepresentatives(ofs, names, assignments, sep);
	}
	else {
		n_total_clusters = graph.saveAssignments(ofs, names, assignments, sep);
	}

	auto dt = std::chrono::high_resolution_clock::now() - t;
	LOG_NORMAL << endl
		<< "  total clusters (including singletons): " << n_total_clusters << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
}