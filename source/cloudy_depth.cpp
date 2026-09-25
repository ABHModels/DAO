#include "cddefines.h"
#include "cddrive.h"
#include "cloudy_depth.h"
#include "compton_cross_section.h"
#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <fcntl.h>
#include <type_traits>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <thread>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>

#undef fopen

namespace {
// Native binary IPC, used only between children and their parent executable.
// LineRec contains strings: serialize each field, never its object representation.
template<class Stream, class T> void scalar(Stream& stream, T& value)
{
    if constexpr (std::is_same_v<Stream,std::ofstream>)
        stream.write(reinterpret_cast<const char*>(&value),sizeof(value));
    else
        stream.read(reinterpret_cast<char*>(&value),sizeof(value));
}
template<class Stream> void string_value(Stream& stream, std::string& value)
{
    size_t n = value.size(); scalar(stream,n);
    if (n > 1024*1024) throw std::runtime_error("Cloudy worker: invalid string size");
    if constexpr (std::is_same_v<Stream,std::ofstream>) stream.write(value.data(),n);
    else { value.resize(n); stream.read(value.data(),n); }
}
template<class Stream> void cell_result(Stream& stream, int id, RadField& rad,
                                        const RTGrids& g, std::vector<LineRec>& lines)
{
    int cell=id, energies=g.NE;
    scalar(stream,cell); scalar(stream,energies);
    if (cell!=id || energies!=g.NE) throw std::runtime_error("Cloudy worker: wrong result cell/grid");
    for (double* field : {rad.T_K,rad.n_e,rad.heating,rad.cooling,rad.line_heat})
        scalar(stream,field[id]);
    for (double** field : {rad.jnu,rad.kabs,rad.ksct,rad.jnu_line})
        for (int e=0;e<g.NE;++e) scalar(stream,field[id][e]);
    size_t n=lines.size(); scalar(stream,n);
    if (n>10000000) throw std::runtime_error("Cloudy worker: invalid line count");
    if constexpr (!std::is_same_v<Stream,std::ofstream>) lines.resize(n);
    for (auto& line:lines) {
        scalar(stream,line.bin); scalar(stream,line.ip);
        scalar(stream,line.emiss); scalar(stream,line.dE_eV); scalar(stream,line.E_eV);
        scalar(stream,line.wave); scalar(stream,line.kappaL); scalar(stream,line.kappaLo);
        scalar(stream,line.damp); scalar(stream,line.y); scalar(stream,line.redis);
        string_value(stream,line.label); string_value(stream,line.comment);
    }
}

struct Batch {
    std::string dir;
    std::map<pid_t,int> active;
    std::vector<std::string> sed_files;
    explicit Batch(int count) {
        char name[]="/tmp/dao-cloudy-XXXXXX";
        if (!mkdtemp(name)) throw std::runtime_error("Cannot create Cloudy worker directory");
        dir=name;
        for(int id=0;id<count;++id)
            sed_files.push_back(std::filesystem::path(dir).filename().string()+"-cell"+std::to_string(id)+".sed");
    }
    ~Batch() {
        // A failed worker must not leave sibling calculations running.
        for (auto [pid,id]:active) kill(pid,SIGKILL);
        for (auto [pid,id]:active) { while (waitpid(pid,nullptr,0)<0 && errno==EINTR) {} }
        std::error_code ignored;
        for(const auto& file:sed_files) std::filesystem::remove(file,ignored);
        std::filesystem::remove_all(dir,ignored);
    }
    std::string path(int id,const char* suffix) const {
        return dir+"/cell"+std::to_string(id)+suffix;
    }
};

void run_cell(int id, RadField& rad, const RTGrids& g, const ModelParams& par,
              int iteration, std::vector<LineRec>& lines, const std::string& sed)
{
    CloudyInput input(g); input.init(par);
    cdInit(); cdTalk(false);
    input.issue_constant();
    // Cloudy requires table SED files to be in the current directory.
    const std::string local_sed=std::filesystem::path(sed).filename().string();
    input.issue_depth(id,rad,g,par,local_sed.c_str());
    if (cdDrive()) throw std::runtime_error("Cloudy cdDrive failed at depth "+std::to_string(id+1));
    lines.clear();
    extract_cloudy_output(id,rad,g,iteration,lines);
    compute_compton_opacity(rad.ksct[id],g.NE,g.ene,rad.T_K[id],pow(10.0,par.nh));
}

unsigned worker_count(unsigned requested,int count)
{
    if (!requested) {
        requested=std::min(4u,std::max(1u,std::thread::hardware_concurrency()));
        if (const char* value=std::getenv("DAO_CLOUDY_WORKERS")) {
            char* end=nullptr; long n=std::strtol(value,&end,10);
            if (end==value || *end || n<1 || n>64)
                throw std::invalid_argument("DAO_CLOUDY_WORKERS must be 1..64");
            requested=unsigned(n);
        }
    }
    if (requested>64) throw std::invalid_argument("Cloudy worker count must be 1..64");
    return std::min(requested,unsigned(std::max(0,count)));
}
}

void run_cloudy_depths(RadField& rad, const RTGrids& g, const ModelParams& par,
                      int iteration, std::vector<std::vector<LineRec>>& lines,
                      unsigned workers)
{
    workers=worker_count(workers,g.ND_MID);
    lines.resize(g.ND_MID);
    Batch batch(g.ND_MID);
    const auto start=std::chrono::steady_clock::now();
    fprintf(stdout,"  [Cloudy] %d depth cells, %u worker process(es)\n",g.ND_MID,workers);
    fflush(stdout);
    if (workers<=1) {
        for (int id=0;id<g.ND_MID;++id)
            run_cell(id,rad,g,par,iteration,lines[id],batch.sed_files[id]);
    } else {
        int next=0, completed=0;
        while (completed<g.ND_MID) {
            while (next<g.ND_MID && batch.active.size()<workers) {
                const int id=next++;
                // The parent is single-threaded here; RT/kernel threads have joined.
                fflush(nullptr);
                const pid_t pid=fork();
                if (pid<0) throw std::runtime_error("Cannot fork Cloudy worker");
                if (pid==0) {
                    const std::string log=batch.path(id,".log");
                    int fd=open(log.c_str(),O_WRONLY|O_CREAT|O_TRUNC,0600);
                    if(fd<0) _exit(2);
                    dup2(fd,STDOUT_FILENO); dup2(fd,STDERR_FILENO); close(fd);
                    try {
                        run_cell(id,rad,g,par,iteration,lines[id],batch.sed_files[id]);
                        std::ofstream out;
                        out.exceptions(std::ios::failbit|std::ios::badbit);
                        out.open(batch.path(id,".bin"),std::ios::binary);
                        cell_result(out,id,rad,g,lines[id]); out.close();
                        fflush(nullptr); _exit(0);
                    } catch (const std::exception& e) {
                        fprintf(stderr,"Cloudy worker failed: %s\n",e.what());
                    } catch (...) {
                        fprintf(stderr,"Cloudy worker raised an exception\n");
                    }
                    fflush(nullptr); _exit(1);
                }
                batch.active.emplace(pid,id);
            }
            bool finished=false;
            for (auto it=batch.active.begin();it!=batch.active.end();) {
                int status=0;
                const pid_t done=waitpid(it->first,&status,WNOHANG);
                if(done<0) {
                    if(errno==EINTR) { ++it; continue; }
                    throw std::runtime_error("Cannot wait for Cloudy worker");
                }
                if(!done) { ++it; continue; }
                const int id=it->second;
                it=batch.active.erase(it); finished=true;
                if (!WIFEXITED(status) || WEXITSTATUS(status)!=0) {
                    std::ifstream log(batch.path(id,".log"));
                    std::string tail((std::istreambuf_iterator<char>(log)),{});
                    if(tail.size()>16000) tail=tail.substr(tail.size()-16000);
                    throw std::runtime_error("Cloudy worker failed at depth "+std::to_string(id+1)+"\n"+tail);
                }
                ++completed;
                fprintf(stdout,"  [Cloudy] completed %d/%d (depth %d)\n",completed,g.ND_MID,id+1);
                fflush(stdout);
            }
            if(!finished) usleep(10000);
        }
        // Deterministic gather; no partial result reaches line escape or RT.
        for (int id=0;id<g.ND_MID;++id) {
            std::ifstream in;
            in.exceptions(std::ios::failbit|std::ios::badbit);
            in.open(batch.path(id,".bin"),std::ios::binary);
            cell_result(in,id,rad,g,lines[id]);
        }
    }
    for(int id=0;id<g.ND_MID;++id)
        fprintf(stdout,"  depth %3d/%d  tau_ref=%.3e  logT=%.3f  log(I)=%.3f  ne/nh=%.3f  H/C=%.3f\n",
                id+1,g.ND_MID,g.tau_mid[id],log10(rad.T_K[id]),rad.log_xi[id]+par.nh,
                rad.n_e[id]/pow(10.0,par.nh),rad.heating[id]/rad.cooling[id]);
    fprintf(stdout,"  [Cloudy] wall time %.3fs\n",
            std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
}
