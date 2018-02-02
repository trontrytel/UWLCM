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
    
  //"--outfreq=10 --nt=200 --spinup=100 --nx=33 --nz=76 --X=3325 --Z=1495 --dt=1 --relax_th_rv=false";
  //"--outfreq=200 --nt=7200 --spinup=6000 --nx=97 --nz=301 --dt=1 --relax_th_rv=false";
  //"--micro=blk_1m --outdir=out_blk_1m --adv_serial=false --async=true --backend=OpenMP --case=dycoms \
  //"--micro=blk_1m --outdir=out_blk_1m --adv_serial=false --async=true --backend=serial --case=dycoms \

  string opts_common = 
    "--outfreq=100 --nt=2000 --spinup=6000 --nx=97 --nz=301 --dt=1 --relax_th_rv=false";
  set<string> opts_micro({
    "--micro=blk_2m --outdir=out_blk_2m --adv_serial=false --async=true --backend=serial --case=dycoms \
     --slice=true --piggy=true \
     --vel_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/dycoms_velocity.h5' \
     --init_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/dycoms_init.h5' \
     --w_src=0 --uv_src=0 --rv_src=0 --th_src=0 --subsidence=0 \
     --cond=0 --cevp=0 --revp=0 --conv=0 --accr=0 --sedi=0 "
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
