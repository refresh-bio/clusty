// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "console.h"
#include "log.h"

#include "graph_named.h"
#include "graph_numbered.h"
#include "sparse_matrix.h"
#include "io.h"

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

	Params::Status status = params.parse(argc, argv);

	if (status == Params::Status::ShowVersion) {
		LOG_NORMAL << VERSION;
		return false;
	}
	else {

		LOG_NORMAL << "Clusty" << endl
			<< "  version " << VERSION
#ifdef GIT_COMMIT
			<< "-" << TOSTRING(GIT_COMMIT)
#endif
			<< " (" << DATE << ")" << endl << endl;

		if (status == Params::Status::Incorrect) {
			params.printUsage();
			return false;
		}

		if (params.verbose) {
			Log::getInstance(Log::LEVEL_VERBOSE).enable();
			Log::getInstance(Log::LEVEL_DEBUG).enable();
		}

		return true;
	}
}

// *******************************************************************************************
std::unique_ptr<Graph> Console::loadGraph(const Params& params) {
	
	unique_ptr<Graph> graph;

	if (params.numericIds) {
		if (needDistances(params.algo)) {
			graph = make_unique<GraphNumbered<dist_t>>(params.numThreads);
		}
		else {
			graph = make_unique<GraphNumbered<mini_dist_t>>(params.numThreads);
		}
	}
	else {
		if (needDistances(params.algo)) {
			graph = make_unique<GraphNamed<dist_t>>(params.numThreads);
		}
		else {
			graph = make_unique<GraphNamed<mini_dist_t>>(params.numThreads);
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
	std::vector<std::string_view>& names) {

	names.clear();
	objects.resize(graph.getNumVertices());
	std::iota(objects.begin(), objects.end(), 0);
	ifstream ifs;

	if (!params.objectsFile.empty()) {
		LOG_NORMAL << "Loading objects from " << params.objectsFile << "... ";
		
		ifstream ifs;
		ifs.open(params.objectsFile, ios_base::binary);

		if (ifs) {
			auto t = std::chrono::high_resolution_clock::now();
			string token;
			
			// omit header
			getline(ifs, token);
			auto is_sep = [](char c) {return c == ',' || c == '\t' || c == '\r' || c == '\n'; };
			auto is_newline = [](char c) {return c == '\r' || c == '\n'; };

			InputBuffer buf(128ULL << 20);

			bool continueReading = true;
			while (continueReading) {

				size_t n_wanted = buf.data + buf.size - buf.block_begin;
				ifs.read(buf.block_begin, n_wanted);
				size_t n_read = (ifs) ? n_wanted : ifs.gcount();

				// reset block
				buf.block_end = buf.block_begin + n_read;
				buf.block_begin = buf.data;

				// if previous block finished between /r and /n characters
				while (is_newline(*buf.block_begin)) {
					++buf.block_begin;
				}

				int n_tail = 0;

				// no more data
				if (n_read < n_wanted) {
					continueReading = false;
				}
				else {
					// find last newline
					while (!is_newline(*(buf.block_end - 1))) {
						--buf.block_end;
						++n_tail;
					}
				}

				// process buffer
				char* p = buf.block_begin;
				bool reachedNewline = false;

				while (p != buf.block_end) {
					char* q = find_if(p, buf.block_end, is_sep);
					reachedNewline = is_newline(*q);
					*q = 0;

					// store name
					size_t len = q - p;
					
					if (len > 0) {
						char* dst = namesBuffer.resize_for_additional(len + 1);
						std::copy_n(p, len + 1, dst); // 0 is already there	
						names.emplace_back(dst, len);
					}

					p = q;
					
					// find newline if not there already
					if (!reachedNewline) {
						p = find_if(p, buf.block_end, is_newline);
					}

					p = std::find_if(p, buf.block_end, [](char c) { return c != '\r' && c != '\n' && c != 0; });
				}

				// copy remaining part after consuming all the lines
				if (continueReading && n_tail > 0) {
					memcpy(buf.data, buf.block_end, n_tail);
					buf.block_begin = buf.data + n_tail;
					buf.block_end = nullptr;
				}
				else {
					buf.block_begin = buf.data;
					buf.block_end = nullptr;
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
	const std::vector<std::string_view>& names,
	const std::vector<int>& assignments) {

	LOG_NORMAL << "Saving clusters (representatives = " << std::boolalpha << params.outputRepresentatives << ")... ";
	auto t = std::chrono::high_resolution_clock::now();

	char sep = params.outputCSV ? ',' : '\t';

	ofstream ofs(params.output, ios_base::binary);
	int n_total_clusters = 0;
	
	n_total_clusters = graph.saveAssignments(ofs, names, assignments, sep, params.outputRepresentatives);
	
	auto dt = std::chrono::high_resolution_clock::now() - t;
	LOG_NORMAL << endl
		<< "  total clusters (including singletons): " << n_total_clusters << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
}