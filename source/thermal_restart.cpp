#include "thermal_restart.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
[[noreturn]] void invalid(const std::string& message)
{ throw std::runtime_error("Cannot resume: "+message); }

std::string read_text(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if(!in) invalid("missing "+path.string());
    std::ostringstream text; text<<in.rdbuf();
    if(in.bad()) invalid("read error in "+path.string());
    return text.str();
}

// Read a required scalar from DAO's own flat params.json, not arbitrary JSON.
std::string scalar(const std::string& text, const std::string& key)
{
    const std::regex pattern("\""+key+"\"\\s*:\\s*(\"[^\"\\\\]*\"|true|false|[-+0-9.eE]+)\\s*[,}]");
    std::smatch match;
    if(!std::regex_search(text,match,pattern)) invalid("missing/invalid metadata: "+key);
    return match[1];
}

bool has_key(const std::string& text, const std::string& key)
{ return std::regex_search(text,std::regex("\""+key+"\"\\s*:")); }

void require(const std::string& text, const std::string& key, const std::string& expected)
{
    if(scalar(text,key)!=expected) invalid("metadata mismatch: "+key);
}

bool near(double value, double expected, double tolerance=1e-6)
{ return std::isfinite(value) && std::abs(value-expected)<=tolerance*std::max(std::abs(expected),1e-30); }

void require_number(const std::string& text, const std::string& key, double expected)
{
    const std::string token=scalar(text,key);
    size_t used=0;
    const double value=std::stod(token,&used);
    // params.json has six significant digits; saved grids have seven.
    if(used!=token.size() || !near(value,expected,1e-5)) invalid("metadata mismatch: "+key);
}

std::filesystem::path iteration_file(const ModelParams& p, const char* stem)
{
    char name[96];
    std::snprintf(name,sizeof(name),"%s_iter%03d.dat",stem,p.resume_iteration);
    return std::filesystem::path(p.run_dir)/name;
}

bool row(std::ifstream& in, std::vector<double>& values, size_t columns)
{
    std::string line;
    while(std::getline(in,line)) {
        const auto start=line.find_first_not_of(" \t\r");
        if(start==std::string::npos || line[start]=='#') continue;
        std::istringstream data(line);
        values.clear(); double value;
        while(data>>value) values.push_back(value);
        if(!data.eof() || values.size()!=columns ||
           !std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);}))
            invalid("malformed or non-finite saved row");
        return true;
    }
    if(in.bad()) invalid("saved data read error");
    return false;
}

void header(std::ifstream& in, const ModelParams& p, const char* stem)
{
    const auto path=iteration_file(p,stem);
    in.open(path);
    if(!in) invalid("missing "+path.string());
    std::string first; std::getline(in,first);
    const std::regex iteration("iter="+std::to_string(p.resume_iteration)+"(\\s|$)");
    if(!std::regex_search(first,iteration)) invalid("iteration header mismatch: "+path.string());
    if(std::string(stem)=="moments" && first.find("Cell-volume-average")==std::string::npos)
        invalid("moments use an incompatible spatial convention");
}
}

void validate_thermal_resume(const ModelParams& p)
{
    if(p.resume_iteration<=0 || p.test_rt || p.max_outer_iterations<=p.resume_iteration)
        invalid("invalid iteration range or test mode");
    const std::string metadata=read_text(std::filesystem::path(p.run_dir)/"params.json");
    require(metadata,"hash","\""+std::string(p.run_hash)+"\"");
    require(metadata,"corona","\""+std::string(p.corona)+"\"");
    require(metadata,"test_mode","\"none\"");
    if(p.frac>0) require(metadata,"frac_convention","\"disk_over_corona_flux_v1\"");
    if(has_key(metadata,"illumination_normalization"))
        require(metadata,"illumination_normalization","\"reflionx\"");
    else {
        // The benchmark predates the consolidated normalization metadata.
        require(metadata,"reflionx_norm","true");
        require(metadata,"xillver_norm","false");
    }
    if(has_key(metadata,"transfer_solver"))
        require(metadata,"transfer_solver","\"constant_cell_v1\"");
    require(metadata,"intensity_location","\"cell-volume-average\"");
    require(metadata,"solver","\"rt_energy_balance_v1\"");
    require(metadata,"scattering_density","\"total_electrons_1.21_nH\"");
    require(metadata,"temperature_iteration","\"local_cell_response_v3\"");
    require(metadata,"thermal_update","\"full\"");
    require(metadata,"rt_thermal_balance","true");
    require(metadata,"test_rt","false");
    require(metadata,"angsca",p.angsca ? "true" : "false");
    require(metadata,"log_depth","true");
    require_number(metadata,"depth_cells",RTGrids::ND_MID_DEFAULT);
    if(has_key(metadata,"angular_points"))
        require_number(metadata,"angular_points",RTGrids::NA);
    else {
        // Legacy benchmark params omit NA. Validate the actual saved GL nodes
        // from the emergent header instead of assuming a default angle count.
        std::ifstream emergent; header(emergent,p,"emergent");
        std::string line; std::vector<double> nodes;
        const std::regex angle("I\\(mu=([-+0-9.eE]+)\\)");
        while(std::getline(emergent,line) && !line.empty() && line[0]=='#') {
            std::smatch match;
            if(std::regex_search(line,match,angle)) nodes.push_back(std::stod(match[1]));
        }
        RTGrids angles; angles.init_angle();
        if(nodes.size()!=RTGrids::NA) invalid("legacy angle count mismatch");
        for(int m=0;m<RTGrids::NA;++m)
            if(!near(nodes[m],angles.mu[m],3e-6)) invalid("legacy angle node mismatch");
    }
    if(p.i_incidence>=0) {
        // read_params validates the requested cosine; snap_incidence then
        // replaces p.incidence with its negative GL node before loading T/J.
        const double requested=std::stod(scalar(metadata,"incidence"));
        RTGrids angles; angles.init_angle();
        int index=0;
        for(int m=1;m<RTGrids::NA;++m)
            if(std::abs(angles.mu[m]+std::abs(requested))<std::abs(angles.mu[index]+std::abs(requested))) index=m;
        if(index!=p.i_incidence || !near(p.incidence,angles.mu[index],1e-12))
            invalid("snapped incidence mismatch");
    } else require_number(metadata,"incidence",p.incidence);
    for(const auto& entry:std::vector<std::pair<std::string,double>>{
            {"nh",p.nh},{"zeta",p.zeta},{"frac",p.frac},
            {"Afe",p.Afe},{"Gamma",p.Gamma},{"E_cut",p.E_cut},{"E_lo_cut",p.E_lo_cut},
            {"kT_e",p.kT_e},{"kT_bb",p.kT_bb},{"taup",p.taup},{"kT_disk",p.kT_disk},
            {"tau_slab",p.tau_slab},{"maxiter",p.maxiter},{"ktype",p.ktype},
            {"E_rt_lo",p.E_rt_lo},{"E_rt_hi",p.E_rt_hi},
            {"column_energy_tolerance",ModelParams::thermal_column_tolerance}})
        require_number(metadata,entry.first,entry.second);

    std::ifstream history(std::filesystem::path(p.run_dir)/"thermal_history.dat");
    if(!history) invalid("missing thermal_history.dat");
    std::vector<double> values;
    int previous=0;
    while(row(history,values,7)) {
        if(values[0]!=previous+1 || values[5]<=0) invalid("incomplete/invalid thermal history");
        ++previous;
    }
    if(previous!=p.resume_iteration)
        invalid("resume iteration must equal the last recorded iteration ("+std::to_string(previous)+")");

    // Never overwrite already completed spectra, even if the history is stale.
    const std::regex product("(emergent|moments|profile)_iter([0-9]+)\\.dat");
    for(const auto& file:std::filesystem::directory_iterator(p.run_dir)) {
        std::smatch match; const auto name=file.path().filename().string();
        if(std::regex_match(name,match,product) && std::stoll(match[2])>p.resume_iteration)
            invalid("newer output already exists: "+name);
    }
    for(const char* stem:{"moments","profile"}) {
        std::ifstream in; header(in,p,stem);
    }
}

void load_thermal_resume(RadField& rad, const RTGrids& g, const ModelParams& p)
{
    // Validate again immediately before use; no writes occur before completion.
    validate_thermal_resume(p);
    std::ifstream profile; header(profile,p,"profile");
    std::vector<double> values;
    std::vector<std::vector<double>> profiles;
    for(int d=0;d<g.ND_MID;++d) {
        if(!row(profile,values,9) || values[0]!=d || !near(values[1],g.tau_mid[d]) ||
           values[2]<=0 || values[3]<0 || values[4]<0 || values[5]<0 || values[7]<0)
            invalid("profile depth/grid/temperature mismatch at cell "+std::to_string(d));
        profiles.push_back(values);
    }
    if(row(profile,values,9)) invalid("profile has extra cells");

    std::ifstream moments; header(moments,p,"moments");
    std::vector<double> J0(size_t(g.ND_MID)*g.NE), J2(J0.size()), J3(J0.size());
    for(int d=0;d<g.ND_MID;++d) for(int e=0;e<g.NE;++e) {
        if(!row(moments,values,7) || !near(values[0],g.ene[e]) || values[1]!=d ||
           !near(values[2],g.tau_mid[d]) || !near(values[3],profiles[d][2]) || values[4]<0)
            invalid("moments energy/depth/T/J mismatch at cell "+std::to_string(d)+", bin "+std::to_string(e));
        const size_t index=size_t(d)*g.NE+e;
        J0[index]=values[4]; J2[index]=values[5]; J3[index]=values[6];
    }
    if(row(moments,values,7)) invalid("moments have extra cells or energy bins");

    for(int d=0;d<g.ND_MID;++d) {
        rad.T_K[d]=profiles[d][2]; rad.n_e[d]=profiles[d][3];
        rad.heating[d]=profiles[d][4]; rad.cooling[d]=profiles[d][5];
        rad.log_xi[d]=profiles[d][6]; rad.line_heat[d]=profiles[d][7];
        for(int e=0;e<g.NE;++e) {
            const size_t index=size_t(d)*g.NE+e;
            rad.J0[d][e]=J0[index]; rad.J2[d][e]=J2[index]; rad.J3[d][e]=J3[index];
        }
    }
}
