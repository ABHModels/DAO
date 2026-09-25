#ifndef DAO_KERNEL_ROW_SPOOL_H
#define DAO_KERNEL_ROW_SPOOL_H

#include <cstdio>
#include <stdexcept>
#include <limits>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <exception>
#include "kernel_payload.h"

// Retain already evaluated, trimmed kernel rows without holding an additional
// full kernel in RAM. tmpfile removes the scratch file when closed. Values are
// written/read as unchanged doubles; no arithmetic or reduction is reordered.
class KernelRowSpool {
	FILE* file_;
	size_t count_=0;
public:
	KernelRowSpool():file_(std::tmpfile()) {
		if(!file_) throw std::runtime_error("kernel: cannot create temporary row spool");
	}
	~KernelRowSpool() { std::fclose(file_); }
	KernelRowSpool(const KernelRowSpool&)=delete;
	KernelRowSpool& operator=(const KernelRowSpool&)=delete;
	void append(const double* row,size_t count) {
		if(count>std::numeric_limits<size_t>::max()-count_ ||
		   std::fwrite(row,sizeof(double),count,file_)!=count)
			throw std::runtime_error("kernel: temporary row spool write failed");
		count_+=count;
	}
	void read(double* data,size_t count) {
		if(count!=count_ || std::fflush(file_)!=0 || std::fseek(file_,0,SEEK_SET)!=0)
			throw std::runtime_error("kernel: temporary row spool size/seek failed");
		if(std::fread(data,sizeof(double),count,file_)!=count)
			throw std::runtime_error("kernel: temporary row spool read failed");
	}
	// A fresh private view discards the previous temperature's normalized
	// pages after they have been streamed to output. Scratch bytes stay exact.
	double* view(KernelPayload& payload,size_t count) {
		if(count!=count_ || std::fflush(file_)!=0 || ::fseeko(file_,0,SEEK_SET)!=0)
			throw std::runtime_error("kernel: cannot prepare spool mapping");
		if(!count) return nullptr;
		double* result=payload.map(file_,count);
		if(!result) throw std::runtime_error("kernel: cannot map retained rows");
		return result;
	}
	// Compute independent rows concurrently but write in the original row
	// order. Each worker owns one scratch row; no floating-point sum changes.
	template<class Compute>
	void evaluate(size_t rows,size_t energies,Compute compute) {
		unsigned workers=std::min(8u,std::max(1u,std::thread::hardware_concurrency()));
		if(const char* value=std::getenv("DAO_KERNEL_THREADS")) {
			char* end=nullptr; const long requested=std::strtol(value,&end,10);
			if(!end || *end || requested<1 || requested>256)
				throw std::invalid_argument("DAO_KERNEL_THREADS must be 1..256");
			workers=unsigned(requested);
		}
		workers=unsigned(std::min(size_t(workers),rows));
		std::fprintf(stdout,"  Kernel construction: %u workers, ordered row output\n",workers);
		std::atomic<size_t> next{0};
		size_t written=0;
		std::mutex mutex;
		std::condition_variable ready;
		std::exception_ptr failure;
		auto work=[&] {
			try {
				std::vector<double> row(energies);
				for(;;) {
					const size_t index=next.fetch_add(1);
					if(index>=rows) break;
					int lo=0,hi=-1;
					compute(index,row.data(),lo,hi);
					std::unique_lock<std::mutex> lock(mutex);
					ready.wait(lock,[&]{return failure || written==index;});
					if(failure) break;
					if(hi>=lo) append(row.data()+lo,size_t(hi-lo+1));
					++written;
					lock.unlock(); ready.notify_all();
				}
			} catch(...) {
				std::lock_guard<std::mutex> lock(mutex);
				if(!failure) failure=std::current_exception();
				ready.notify_all();
			}
		};
		std::vector<std::thread> threads;
		try { for(unsigned i=1;i<workers;++i) threads.emplace_back(work); }
		catch(...) {
			{ std::lock_guard<std::mutex> lock(mutex); failure=std::current_exception(); }
			ready.notify_all();
			for(auto& thread:threads) thread.join();
			std::rethrow_exception(failure);
		}
		work();
		for(auto& thread:threads) thread.join();
		if(failure) std::rethrow_exception(failure);
	}
};
#endif
