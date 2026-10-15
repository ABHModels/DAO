#ifndef DAO_RT_PARALLEL_H
#define DAO_RT_PARALLEL_H

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// Only independent output cells may be distributed here. Each cell keeps its
// original serial integration order; join is a barrier before the next stage.
template<class Work>
void rt_parallel_depths(int count, Work work,
                        const char* thread_variable="DAO_RT_THREADS",
                        unsigned default_limit=16)
{
	unsigned workers=std::min(default_limit,std::max(1u,std::thread::hardware_concurrency()));
	if(const char* value=std::getenv(thread_variable)) {
		char* end=nullptr;
		const long requested=std::strtol(value,&end,10);
		if(!end || *end || requested<1 || requested>256)
			throw std::invalid_argument(std::string(thread_variable)+" must be 1..256");
		workers=unsigned(requested);
	}
	workers=std::min(workers,unsigned(std::max(0,count)));
	if(workers<=1) {
		for(int i=0;i<count;++i) work(i);
		return;
	}
	std::atomic<int> next{0};
	std::atomic<bool> stop{false};
	std::mutex mutex;
	std::exception_ptr failure;
	auto run=[&] {
		try {
			while(!stop.load(std::memory_order_relaxed)) {
				const int i=next.fetch_add(1,std::memory_order_relaxed);
				if(i>=count) break;
				work(i);
			}
		} catch(...) {
			std::lock_guard<std::mutex> lock(mutex);
			if(!failure) failure=std::current_exception();
			stop.store(true,std::memory_order_relaxed);
		}
	};
	std::vector<std::thread> threads;
	try {
		for(unsigned i=1;i<workers;++i) threads.emplace_back(run);
	} catch(...) {
		stop.store(true,std::memory_order_relaxed);
		for(auto& thread:threads) thread.join();
		throw;
	}
	run();
	for(auto& thread:threads) thread.join();
	if(failure) std::rethrow_exception(failure);
}
#endif
