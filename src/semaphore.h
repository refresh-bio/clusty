#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include <mutex>
#include <condition_variable>

// *****************************************************************************************
//
class Semaphore {
protected:
	int counter;
	std::mutex mutex;
	std::condition_variable cv;

public:
	// *****************************************************************************************
	//
	Semaphore() : counter(0) {}

	// *****************************************************************************************
	//
	void inc() {
		std::unique_lock<std::mutex> lk(mutex);
		++counter;
	}

	// *****************************************************************************************
	//
	void inc(int num) {
		std::unique_lock<std::mutex> lk(mutex);
		counter += num;
	}

	// *****************************************************************************************
	//
	void dec() {
		std::unique_lock<std::mutex> lk(mutex);
		--counter;

		if (counter == 0)
			cv.notify_one();
	}

	// *****************************************************************************************
	//
	void dec_notify_all() {
		std::unique_lock<std::mutex> lk(mutex);
		--counter;

		if (counter == 0)
			cv.notify_all();
	}

	// *****************************************************************************************
	//
	void waitForZero() {
		std::unique_lock<std::mutex> lk(mutex);
		cv.wait(lk, [this] {return counter == 0; });
	}
};
