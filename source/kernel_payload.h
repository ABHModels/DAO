#ifndef DAO_KERNEL_PAYLOAD_H
#define DAO_KERNEL_PAYLOAD_H
#include <cstdio>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// Private file mapping preserves the existing writable data API (e.g. the
// standalone normalizer) without ever writing changes into the input cache.
// Normal RT reads fault in only the pages it uses; no full-size heap copy.
class KernelPayload {
	void* base_=MAP_FAILED;
	size_t length_=0;
public:
	KernelPayload()=default;
	KernelPayload(const KernelPayload&)=delete;
	KernelPayload& operator=(const KernelPayload&)=delete;
	~KernelPayload() { unmap(); }
	void unmap() {
		if(base_!=MAP_FAILED) ::munmap(base_,length_);
		base_=MAP_FAILED; length_=0;
	}
	double* map(FILE* file,size_t count) {
		if(!count || count>std::numeric_limits<size_t>::max()/sizeof(double)) return nullptr;
		const off_t offset=::ftello(file);
		struct stat info;
		const long page=::sysconf(_SC_PAGESIZE);
		const size_t bytes=count*sizeof(double);
		if(offset<0 || page<=0 || ::fstat(::fileno(file),&info)!=0 ||
		   offset>info.st_size || bytes>static_cast<unsigned long long>(info.st_size-offset)) return nullptr;
		const off_t aligned=offset-offset%page;
		const size_t skip=static_cast<size_t>(offset-aligned);
		if(skip>std::numeric_limits<size_t>::max()-bytes || offset%alignof(double)!=0) return nullptr;
		void* mapping=::mmap(nullptr,bytes+skip,PROT_READ|PROT_WRITE,MAP_PRIVATE,::fileno(file),aligned);
		if(mapping==MAP_FAILED) return nullptr;
		unmap(); base_=mapping; length_=bytes+skip;
		return reinterpret_cast<double*>(static_cast<char*>(mapping)+skip);
	}
	void release(double*& data) {
		if(base_!=MAP_FAILED) unmap(); else delete[] data;
		data=nullptr;
	}
};

// Never truncate an inode that may still be mapped by this or another run.
// Write in the destination directory and publish only a complete cache.
class KernelCacheOutput {
	std::string target_;
	std::vector<char> temporary_;
	FILE* file_=nullptr;
public:
	explicit KernelCacheOutput(const char* target):target_(target) {
		const std::string pattern=target_+".tmp.XXXXXX";
		temporary_.assign(pattern.begin(),pattern.end()); temporary_.push_back('\0');
		const int fd=::mkstemp(temporary_.data());
		if(fd<0) throw std::runtime_error("kernel: cannot create cache output");
		file_=::fdopen(fd,"w+b");
		if(!file_) { ::close(fd); ::unlink(temporary_.data()); throw std::runtime_error("kernel: cannot open cache output"); }
	}
	~KernelCacheOutput() {
		if(file_) std::fclose(file_);
		if(!temporary_.empty()) ::unlink(temporary_.data());
	}
	KernelCacheOutput(const KernelCacheOutput&)=delete;
	KernelCacheOutput& operator=(const KernelCacheOutput&)=delete;
	FILE* file() const { return file_; }
	double* map_payload(KernelPayload& payload,size_t count,off_t offset) {
		if(std::ferror(file_) || std::fflush(file_)!=0 || ::fseeko(file_,offset,SEEK_SET)!=0)
			throw std::runtime_error("kernel: cannot flush completed payload");
		if(!count) return nullptr;
		double* result=payload.map(file_,count);
		if(!result) throw std::runtime_error("kernel: cannot map completed payload");
		return result;
	}
	void commit() {
		const bool failed=std::ferror(file_)!=0;
		const int result=std::fclose(file_); file_=nullptr;
		if(failed || result!=0 || ::rename(temporary_.data(),target_.c_str())!=0)
			throw std::runtime_error("kernel: failed to publish complete cache");
		temporary_.clear();
	}
};
#endif
