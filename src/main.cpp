// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <fstream>
#include <vector>
#include <memory>
#include <numeric>
#include <iostream>
#include <chrono>
#include <stdexcept>

#include "linkage_heaptrix.h"
#include "uclust.h"
#include "set_cover.h"
#include "single_bfs.h"
#include "cd_hit.h"
#include "leiden.h"

#include "console.h"
#include "log.h"

#define VAL(str) #str
#define TOSTRING(str) VAL(str)
#include "version.h"

using namespace std;

#define debug(x) std::cerr << __FILE__ << " (" << __LINE__ << ") " << #x << " == " << (x) << std::endl


int main(int argc, char** argv) 
{
   
	Log::getInstance(Log::LEVEL_NORMAL).enable();

	try {
	   LOG_NORMAL << "Clusty" << endl
		   << "  version " << VERSION
#ifdef GIT_COMMIT
		   << "-" << TOSTRING(GIT_COMMIT)
#endif
		   << " (" << DATE << ")" << endl << endl;

	Console console;

	if (!console.parse(argc, argv)) {
		console.printUsage();
		return 0;
	}

	shared_ptr<SparseMatrix> dists;
	std::shared_ptr<IClustering<decltype(*dists)>> clustering;

	switch (console.algo)
	{
	case Algo::SingleLinkage:
		clustering = make_shared<SingleLinkageBFS<decltype(*dists)>>(); break;
	case Algo::CompleteLinkage:
		clustering = make_shared<CompleteLinkage<decltype(*dists)>>(); break;
	case Algo::UClust:
		clustering = make_shared<UClust<decltype(*dists)>>(); break;
	case Algo::SetCover:
		clustering = make_shared<SetCover<decltype(*dists)>>(); break;
	case Algo::CdHit:
		clustering = make_shared<CdHit<decltype(*dists)>>(); break;
	case Algo::Leiden:
		clustering = make_shared<Leiden<decltype(*dists)>>(console.leidenParams); break;
	default:
		throw std::runtime_error("Unkown clustering algorithm"); 
	}

	LOG_NORMAL << "Loading pairwise distances from " << console.distancesFile << "... ";
	auto t = std::chrono::high_resolution_clock::now();
	
	vector<char> filebuf(128ULL << 20);  // 128MB buffer
	ifstream ifs;
	ifs.rdbuf()->pubsetbuf(filebuf.data(), filebuf.size());
	ifs.open(console.distancesFile, ios_base::binary);
	
	if (!ifs) {
		throw std::runtime_error("Unable to open distance file");
	}

	if (console.numericIds) {
		dists = make_shared<SparseMatrixNumbered>();
	}
	else {
		dists = make_shared<SparseMatrixNamed>();
	}

	map<DistanceSpecification, distance_transformation_t> transforms {
		{ DistanceSpecification::Distance, [](double d) { return d; } },
		{ DistanceSpecification::Similarity,		[](double d) { return 1.0 - d; } },
		{ DistanceSpecification::PercentSimilarity, [](double d) { return 1.0 - d * 0.01; } },
	};

	size_t n_total_dists = dists->load(ifs, console.idColumns, console.distanceColumn,
		transforms[console.distanceSpecification], console.columns2filters);
	auto dt = std::chrono::high_resolution_clock::now() - t;

	ifs.close();
	LOG_NORMAL << endl
		<< "  input graph: " << dists->num_input_objects() << " nodes, " << n_total_dists << " edges" << endl
		<< "  filtered graph: " << dists->num_objects() << " nodes, " << dists->num_elements() << " edges" << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;

	vector<string> names;
	vector<int> objects(dists->num_objects());
	std::iota(objects.begin(), objects.end(), 0);

	if (!console.objectsFile.empty()) {
		LOG_NORMAL << "Loading objects from " << console.objectsFile << "... ";
		ifs.open(console.objectsFile);
		
		if (ifs) {
			t = std::chrono::high_resolution_clock::now();
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

			// reorder sequences according to the objects file 
			if (auto m = dynamic_pointer_cast<SparseMatrixNamed>(dists)) {
				int obj_id = 0;
				for (const string& name : names) {
					int local_id = m->get_id(name.c_str());
					if (local_id != -1) {
						objects[obj_id++] = local_id;
					}
				}
			}
			else {
				auto q = dynamic_pointer_cast<SparseMatrixNumbered>(dists);
				int obj_id = 0;
				for (size_t i = 0; i < names.size(); ++i) {
					int local_id = q->get_local_id(i);
					if (local_id != -1) {
						objects[obj_id++] = local_id;
					}
				}
			}

			dt = std::chrono::high_resolution_clock::now() - t;
			LOG_NORMAL << "  total objects: " << names.size() << endl
				<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
		}
		else {
			throw std::runtime_error("Unable to open objects file");
		}
	}
	
	LOG_NORMAL << "Clustering... ";
	vector<int> assignments;
	t = std::chrono::high_resolution_clock::now();

	double threshold = std::nexttoward(std::numeric_limits<double>::max(), 0.0);
	int n_clusters = (*clustering)(*dists, objects, threshold, assignments);
	dt = std::chrono::high_resolution_clock::now() - t;
	LOG_NORMAL << endl
		<< "  objects: " << dists->num_objects() << ", clusters: " << n_clusters << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
	
	LOG_NORMAL << "Saving clusters... ";
	t = std::chrono::high_resolution_clock::now();

	char sep = console.outputCSV ? ',' : '\t';

	ofstream ofs(console.output);
	int n_total_clusters = 0;
	if (console.outputRepresentatives) {
		n_total_clusters = dists->saveRepresentatives(ofs, names, assignments, sep);
	}
	else {
		n_total_clusters = dists->saveAssignments(ofs, names, assignments, sep);
	}

	dt = std::chrono::high_resolution_clock::now() - t;
	LOG_NORMAL << endl
		<< "  total clusters (including singletons): " << n_total_clusters << endl
		<< "  time [s]: " << chrono::duration<double>(dt).count() << endl;
   }
   catch  (const std::string & message)
   {
	  LOG_NORMAL << message << std::endl;
	  return -1;
   }
   catch (const std::exception & exc)
   {
	   LOG_NORMAL << exc.what() << std::endl;
	  return -1;
   }
   catch (...)
   {
	  LOG_NORMAL << "Something went wrong and it should not have gone.";
	  return -1;
   }

	return 0;
}
