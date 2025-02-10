#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "graph.h"
#include "sparse_matrix.h"
#include "log.h"
#include "parallel-queues.h"
#include "io.h"
#include "semaphore.h"


#include <thread>
#include <vector>
#include <atomic>
#include <algorithm>
#include <barrier>
#include <thread>

// *******************************************************************************************/
template <class Distance>
class GraphSparse : public Graph {

protected:
	SparseMatrix<Distance> matrix;

public:

	GraphSparse(int numThreads) : Graph(numThreads) {}

	IMatrix& getMatrix() override { return matrix; }

	size_t getNumVertices() const override { return matrix.num_objects(); }

	size_t getNumEdges() const override { return matrix.num_elements(); }

	size_t load(
		std::ifstream& ifs,
		const std::pair<std::string, std::string>& idColumns,
		const std::string& distanceColumn,
		distance_transformation_t transform,
		const std::map<std::string, ColumnFilter>& columns2filters) override;

protected:

	virtual IEdgesCollection* createEdgesCollection(size_t preallocSize) = 0;

	virtual void initLoad();

	virtual void finalizeLoad();

	virtual bool parseBlock(
		char* block_begin,
		char* block_end,
		distance_transformation_t transform,
		IEdgesCollection& edges,
		size_t& n_rows) = 0;

	virtual void updateMappings(
		IEdgesCollection& edges) = 0;

	virtual void extendMatrix() = 0;

	virtual void updateMatrix(
		const IEdgesCollection& edges,
		int startRow,
		int stride) = 0;

};


/*********************************************************************************************************************/
template <class Distance>
void GraphSparse<Distance>::initLoad() {

	// assume space for 8M objects
	matrix.distances.reserve(8LL << 20);
}


/*********************************************************************************************************************/
template <class Distance>
void GraphSparse<Distance>::finalizeLoad() {

	// if neccessary, sort distances in rows according to the second id
	matrix.n_elements = 0;

	std::vector<std::thread> workers(this->numThreads);
	std::atomic<size_t> total_elements = 0;

	for (int tid = 0; tid < numThreads; ++tid) {
		workers[tid] = std::thread([this, tid, &total_elements]() {
			size_t local_elements = 0;
			int n_rows = (int)this->matrix.distances.size();

			for (int i = tid; i < n_rows; i += this->numThreads) {
				auto& row = this->matrix.distances[i];
				std::sort(row.begin(), row.end());
				auto newEnd = std::unique(row.begin(), row.end(), [](const Distance& a, const Distance& b) { return a.get_id() == b.get_id(); });

				row.erase(newEnd, row.end());
				local_elements += row.size();
			}

			total_elements += local_elements;
			});
	}

	for (auto& w : workers) {
		w.join();
	}

	matrix.n_elements = total_elements;

	// Print distance histogram in the verbose mode
	if (Log::getInstance(Log::LEVEL_VERBOSE).isEnabled()) {

		std::vector<double> histo_bounds{ 0 };
		double width = 0.001;

		while (histo_bounds.back() < 0.05)
		{
			histo_bounds.push_back(histo_bounds.back() + width);
		}
		histo_bounds.push_back(std::numeric_limits<double>::max());
		std::vector<int> histo(histo_bounds.size());

		for (auto& row : matrix.distances) {
			for (const auto& e : row) {
				for (size_t i = 0; i < histo_bounds.size(); ++i) {
					if (e.get_d() < histo_bounds[i]) {
						++histo[i];
						break;
					}
				}
			}
		}

		LOG_VERBOSE << std::endl << "Distance histogram" << std::endl;
		for (size_t i = 0; i < histo_bounds.size(); ++i) {
			LOG_VERBOSE << "  d < " << histo_bounds[i] << ": " << histo[i] << std::endl;
		}
		LOG_VERBOSE << std::endl;
	}
}


/*********************************************************************************************************************/
template <class Distance>
size_t GraphSparse<Distance>::load(
	std::ifstream& ifs,
	const std::pair<std::string, std::string>& idColumns,
	const std::string& distanceColumn,
	distance_transformation_t transform,
	const std::map<std::string, ColumnFilter>& columns2filters) {

	std::atomic<size_t> n_total_distances = 0;

	int numParsers = std::max(1, (numThreads - 2) / 2); // at least one parser
	int numUpdaters = std::max(1, (numThreads - 2) / 2); // at least one updater

	// create a vector of input buffers and edges collections
	std::vector<InputBuffer*> buffers(numParsers + 2);
	std::vector<IEdgesCollection*> edgesCollections(numParsers + 2);

	for (int i = 0; i < (int)buffers.size(); ++i) {
		buffers[i] = new InputBuffer(128ULL << 20);
		edgesCollections[i] = this->createEdgesCollection(1ULL << 20);
	}

	struct task_t {
		int buffer_id{ -1 };
		bool buffer_released{ false };
		int collection_id{ -1 };
		int portion_id{ -1 };
	};

	// create queues
	refresh::parallel_queue<int> freeBuffersQueue(buffers.size(), numParsers, "free-buffers-queue");
	refresh::parallel_queue<task_t> blocksQueue(buffers.size(), 1, "blocks-queue");

	refresh::parallel_queue<int> freeCollectionsQueue(edgesCollections.size(), numUpdaters, "free-collections-queue");
	refresh::parallel_priority_queue<task_t> edgesQueue(edgesCollections.size(), numParsers, "edges-queue");

	std::vector<refresh::parallel_queue<task_t>*> updatersQueues;
	for (int i = 0; i < numUpdaters; ++i) {
		updatersQueues.push_back(new refresh::parallel_queue<task_t>(1, 1, "updater-queue-" + std::to_string(i)));
	}

	// get header
	this->processHeader(ifs, idColumns, distanceColumn, columns2filters);

	this->initLoad();

	// add free buffers and edges collecion to queue
	for (int i = 0; i < (int)buffers.size(); ++i) {
		freeBuffersQueue.push(int{ i });
		freeCollectionsQueue.push(int{ i });
	}

	// start parsers
	std::vector<std::thread> parsers(numParsers);
	for (int tid = 0; tid < (int)parsers.size(); ++tid) {
		parsers[tid] = std::thread([tid, this,
			&buffers, &freeBuffersQueue, &blocksQueue, &freeCollectionsQueue, &edgesCollections, &edgesQueue,
			transform, &n_total_distances] () {

				int collection_id;
				task_t task;

				while (freeCollectionsQueue.pop(collection_id) && blocksQueue.pop(task)) {
					LOG_DEBUG << "parser-" << tid << " pop " << task.portion_id << " (buf " << task.buffer_id << ")" << std::endl;

					InputBuffer* buf = buffers[task.buffer_id];
					IEdgesCollection* edges = edgesCollections[collection_id];

					edges->clear();
					size_t n_local_rows = 0;

					bool can_release = this->parseBlock(buf->block_begin, buf->block_end, transform, *edges, n_local_rows);
					n_total_distances += n_local_rows;

					// fill some info
					task.buffer_released = can_release;
					task.collection_id = collection_id;

					LOG_DEBUG << "parser-" << tid << " push " << task.portion_id << "[col " << task.collection_id << "]" << std::endl;
					edgesQueue.push(task.portion_id, task_t{ task });

					if (can_release) {
						LOG_DEBUG << "parser-" << tid << " free(buf " << task.buffer_id << ")" << std::endl;
						freeBuffersQueue.push(int{ task.buffer_id });
					}
				}

				edgesQueue.mark_completed();
			});
	}


	// start mapper
	Semaphore activeUpdaters;
	std::thread mapper([this, &edgesCollections, &freeCollectionsQueue, &edgesQueue, &freeBuffersQueue, &updatersQueues, &activeUpdaters]() {

		task_t task;

		while (edgesQueue.pop(task)) {
			LOG_DEBUG << "mapper pop " << task.portion_id << " [col " << task.collection_id << "]" << std::endl;
			auto edges = edgesCollections[task.collection_id];
			this->updateMappings(*edges);

			// wait with extension until updaters finish previous portion
			activeUpdaters.waitForZero();
			this->extendMatrix();

			// push task to all updaters
			LOG_DEBUG << "mapper push " << task.portion_id << " [col " << task.collection_id << "]" << std::endl;
			activeUpdaters.inc(updatersQueues.size());
			for (auto q : updatersQueues) {
				q->push(task_t{ task });
			}

			if (!task.buffer_released) {
				LOG_DEBUG << "mapper free (buf " << task.buffer_id << ")" << std::endl;
				freeBuffersQueue.push(int{ task.buffer_id });
			}

		}

		for (auto q : updatersQueues) {
			q->mark_completed();
		}
		});

	// start updaters
	std::vector<std::thread> updaters(numUpdaters);
	std::barrier syncPoint(updaters.size());

	for (int tid = 0; tid < updaters.size(); ++tid) {
		updaters[tid] = std::thread([this, tid, &edgesCollections, &freeCollectionsQueue, &updatersQueues, &syncPoint, &activeUpdaters]() {

			task_t task;

			while (updatersQueues[tid]->pop(task)) {

				LOG_DEBUG << "updater-" << tid << " pop " << task.portion_id << "[col " << task.collection_id << "]" << std::endl;
				IEdgesCollection* edges = edgesCollections[task.collection_id];
				this->updateMatrix(*edges, tid, (int)updatersQueues.size());

				// decrement and wait until all updaters finish
				activeUpdaters.dec();
				syncPoint.arrive_and_wait();
				//activeUpdaters.waitForZero();

				if (tid == 0) {
					LOG_DEBUG << "updater-" << tid << " free[col " << task.collection_id << "]" << std::endl;
					freeCollectionsQueue.push(int{ task.collection_id });
				}
			}
			});
	}


	// start loader
	int buffer_id = -1;
	freeBuffersQueue.pop(buffer_id);
	LOG_DEBUG << "loader reserve (buf" << buffer_id << ")" << std::endl;

	bool continueReading = true;
	for (int i_block = 0; continueReading; ++i_block) {

		InputBuffer& buf{ *buffers[buffer_id] };

		size_t n_wanted = buf.data + buf.size - buf.block_begin;
		ifs.read(buf.block_begin, n_wanted);
		size_t n_read = (ifs) ? n_wanted : ifs.gcount();

		// reset block
		buf.block_end = buf.block_begin + n_read;
		buf.block_begin = buf.data;

		int n_tail = 0;

		// no more data
		if (n_read < n_wanted) {
			continueReading = false;
		}
		else {
			// find last newline
			while (!isNewline(*(buf.block_end - 1))) {
				--buf.block_end;
				++n_tail;
			}
		}

		// pop next free buffer
		int next_buffer_id = -1;
		freeBuffersQueue.pop(next_buffer_id);
		InputBuffer& nextBuf{ *buffers[next_buffer_id] };

		LOG_DEBUG << "loader reserve (buf " << next_buffer_id << ")" << std::endl;

		// copy remaining part after consuming all the lines
		if (continueReading && n_tail > 0) {
			memcpy(nextBuf.data, buf.block_end, n_tail);
			nextBuf.block_begin = nextBuf.data + n_tail;
			nextBuf.block_end = nullptr;
		}
		else {
			nextBuf.block_begin = nextBuf.data;
			nextBuf.block_end = nullptr;
		}

		LOG_DEBUG << "loader push " << i_block << " (buf " << buffer_id << ")" << std::endl;
		blocksQueue.push(task_t{ buffer_id, false, -1, i_block });
		buffer_id = next_buffer_id;
	}

	blocksQueue.mark_completed();

	// join threads
	for (auto& t : parsers) { t.join(); }
	mapper.join();
	for (auto& t : updaters) { t.join(); }

	this->finalizeLoad();

	// free memory 
	for (auto& e : buffers) { delete e; }
	for (auto& e : edgesCollections) { delete e; }
	for (auto& e : updatersQueues) { delete e; }

	return n_total_distances;
}