#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"
#include "bins.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string bins_wet_str, bins_dry_str, outdir;

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_dry();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_dry_str = tmp.str();
  }

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_wet();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_wet_str = tmp.str();
  }

  string opts_common = 
    "--outfreq=200 --nt=12000 --spinup=9600 --nx=97 --nz=301 --dt=0.75 --relax_th_rv=false"; // DYCOMS: 128x300 ; 600 21600 3600
  set<string> opts_micro({
    "--adv_serial=false --async=true --micro=lgrngn --outdir=out_lgrngn_s_2017 --backend=CUDA --sd_conc=2048 --sstp_cond=1 --sstp_coal=1 --case=dycoms --rng_seed=2017 "
      " --out_wet=\""
        ".5e-6:25e-6|0,1,2,3,4,5,6;" // FSSP
        "25e-6:1|0,1,2,3,4,5,6;"     // "rain"
        ".5e-6:1|0,1,2,3,4,5,6;"     // all hydro
       + bins_wet_str +  // aerosol spectrum (wet)
        "\""
      " --out_dry=\""
        + bins_dry_str + // aerosol spectrum (dry)
      "\""
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m;  
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}
