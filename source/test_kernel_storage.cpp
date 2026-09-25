#include "kernel_payload.h"
#include "kernel_row_spool.h"
#include "rt_parallel.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>

int main() {
	char directory[]="/tmp/dao-kernel-storage.XXXXXX";
	assert(::mkdtemp(directory));
	const std::string path=std::string(directory)+"/payload.bin";
	const double original[]={1.,0.,-0.,1e-200,1e200};
	{
		KernelCacheOutput out(path.c_str());
		assert(std::fwrite(original,sizeof(double),5,out.file())==5);
		out.commit();
	}
	KernelPayload mapping;
	FILE* fp=std::fopen(path.c_str(),"rb"); assert(fp);
	double* data=mapping.map(fp,5); assert(data);
	std::fclose(fp);
	assert(std::memcmp(data,original,sizeof(original))==0);
	data[0]=2.; // private mapping must not modify the stored input
	double disk[5];
	fp=std::fopen(path.c_str(),"rb"); assert(fp);
	assert(std::fread(disk,sizeof(double),5,fp)==5); std::fclose(fp);
	assert(std::memcmp(disk,original,sizeof(original))==0);
	{
		KernelCacheOutput out(path.c_str());
		assert(std::fwrite(data,sizeof(double),5,out.file())==5);
		out.commit(); // replace the filename while its old inode remains mapped
	}
	assert(data[0]==2. && data[4]==original[4]);
	fp=std::fopen(path.c_str(),"rb"); assert(fp);
	assert(std::fread(disk,sizeof(double),5,fp)==5); std::fclose(fp);
	assert(disk[0]==2. && std::memcmp(disk+1,original+1,4*sizeof(double))==0);
	{
		KernelCacheOutput abandoned(path.c_str());
		std::fputs("incomplete",abandoned.file()); // no commit: preserve target
	}
	fp=std::fopen(path.c_str(),"rb"); assert(fp);
	KernelPayload invalid;
	assert(invalid.map(fp,6)==nullptr); // file too short
	assert(std::fread(disk,sizeof(double),5,fp)==5 && disk[0]==2.);
	std::fclose(fp);
	mapping.release(data); assert(data==nullptr);
	mapping.release(data); // harmless repeated cleanup
	{
		KernelPayload scratch;
		KernelRowSpool spool;
		spool.append(original,5);
		double* values=spool.view(scratch,5);
		assert(std::memcmp(values,original,sizeof(original))==0);
		values[0]=42.;
		values=spool.view(scratch,5); // discard previous private changes
		assert(std::memcmp(values,original,sizeof(original))==0);
		KernelCacheOutput output(path.c_str());
		const double header=17.;
		assert(std::fwrite(&header,sizeof(double),1,output.file())==1);
		assert(std::fwrite(values,sizeof(double),5,output.file())==5);
		values=output.map_payload(scratch,5,sizeof(double));
		output.commit();
		assert(std::memcmp(values,original,sizeof(original))==0);
		scratch.release(values);
	}
	{
		KernelPayload empty;
		KernelRowSpool spool;
		assert(spool.view(empty,0)==nullptr);
		bool failed=false;
		try { spool.view(empty,1); }
		catch(const std::runtime_error&) { failed=true; }
		assert(failed);
	}
	{
		KernelRowSpool spool;
		spool.append(original,2); spool.append(original+2,3);
		spool.read(disk,5);
		assert(std::memcmp(disk,original,sizeof(original))==0);
	}
	::setenv("DAO_KERNEL_THREADS","3",1);
	{
		KernelRowSpool spool;
		spool.evaluate(100,4,[](size_t index,double* row,int& lo,int& hi) {
			lo=1; hi=2; row[1]=double(index); row[2]=-double(index);
		});
		double rows[200]; spool.read(rows,200);
		for(int i=0;i<100;++i) assert(rows[2*i]==i && rows[2*i+1]==-i);
	}
	{
		KernelRowSpool spool;
		bool failed=false;
		try {
			spool.evaluate(100,4,[](size_t index,double* row,int& lo,int& hi) {
				if(index==2) throw std::runtime_error("injected worker failure");
				lo=hi=0; row[0]=double(index);
			});
		} catch(const std::runtime_error&) { failed=true; }
		assert(failed); // all workers must join rather than deadlock
	}
	::unsetenv("DAO_KERNEL_THREADS");
	::setenv("DAO_RT_THREADS","3",1);
	{
		std::vector<int> visits(101,0);
		rt_parallel_depths(int(visits.size()),[&](int i) { ++visits[i]; });
		for(int count:visits) assert(count==1);
		bool failed=false;
		try {
			rt_parallel_depths(101,[](int i) {
				if(i==2) throw std::runtime_error("injected RT worker failure");
			});
		} catch(const std::runtime_error&) { failed=true; }
		assert(failed);
		rt_parallel_depths(0,[](int) { assert(false); });
	}
	::unsetenv("DAO_RT_THREADS");
	assert(::unlink(path.c_str())==0);
	assert(::rmdir(directory)==0); // also verifies abandoned temporary cleanup
	std::puts("PASS: exact spool, private mapping, atomic replacement, failed-output cleanup, truncated payload rejection, RT worker coverage and exception cleanup.");
}
