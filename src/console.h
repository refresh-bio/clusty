// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************
#pragma once
#include "params.h"
#include "graph.h"

#include "linkage_heaptrix.h"
#include "uclust.h"
#include "set_cover.h"
#include "single_bfs.h"
#include "cd_hit.h"
#include "leiden.h"

#include <memory>
#include <vector>
#include <string>

class Console {
	Params params;

public:		
	bool init(int argc, char** argv, Params& params);
	
	std::unique_ptr<Graph> loadGraph(const Params& params);
	
	void loadObjects(
		const Params& params,
		const Graph& graph,
		std::vector<int>& objects,
		std::vector<std::string>& names);

	void doClustering(
		const Params& params,
		Graph& graph,
		const std::vector<int>& objects,
		std::vector<int>& assignments);

	void saveAssignments(
		const Params& params,
		const Graph& graph,
		const std::vector<std::string>& names,
		const std::vector<int>& assignments);

protected:

	bool needDistances(Algo algo) const { 
		return (algo == Algo::CompleteLinkage || algo == Algo::Leiden || algo == Algo::UClust);
	}

	template <class Distance>
	std::unique_ptr<IClustering<Distance>> createClusteringAlgo(const Params& params) {
		
		std::unique_ptr<IClustering<Distance>> clustering;

		switch (params.algo)
		{
			
		case Algo::SingleLinkage:
			clustering = std::make_unique<SingleLinkageBFS<Distance>>(); break;
		case Algo::CompleteLinkage:
			clustering = std::make_unique<CompleteLinkage<Distance>>(); break;
		case Algo::UClust:
			clustering = std::make_unique<UClust<Distance>>(); break;
		case Algo::SetCover:
			clustering = std::make_unique<SetCover<Distance>>(); break;
		case Algo::CdHit:
			clustering = std::make_unique<CdHit<Distance>>(); break;
		case Algo::Leiden:
			clustering = std::make_unique<Leiden<Distance>>(params.leidenParams); break;
			
		default:
			throw std::runtime_error("Unkown clustering algorithm");
		}

		return clustering;
	}

};
