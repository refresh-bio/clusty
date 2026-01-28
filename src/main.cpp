// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#include "console.h"
#include "log.h"

#include <vector>
#include <stdexcept>

using namespace std;

int main(int argc, char** argv)
{
	try {

		Console console;
		Params params;

		std::vector<int> objects;
		std::vector<string_view> names;
		std::vector<int> assignments;

		if (!console.init(argc, argv, params)) {
			return 0;
		}

		std::unique_ptr<Graph> graph = console.loadGraph(params);

		console.loadObjects(params, *graph, objects, names);
		if (graph->getNumEdges() > 0) {
			console.doClustering(params, *graph, objects, assignments);
		}
		else if (graph->getNumVertices() > 0) {
			// No edges but have vertices: each vertex is its own singleton cluster
			assignments.resize(graph->getNumVertices());
			std::iota(assignments.begin(), assignments.end(), 0);
		}
		console.saveAssignments(params, *graph, names, assignments);

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
