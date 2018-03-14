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
    
  //"--outfreq=10  --nt=200  --spinup=100  --nx=33 --nz=76  --X=3325 --Z=1495 --dt=1 --relax_th_rv=false";
  //"--outfreq=200 --nt=9000 --spinup=7200 --nx=97 --nz=301 --dt=1                   --relax_th_rv=false";

  //"--micro=blk_1m --outdir=out_blk_1m --adv_serial=false --async=true --backend=OpenMP --case=dycoms \
  //"--micro=blk_1m --outdir=out_blk_1m --adv_serial=false --async=true --backend=serial --case=dycoms \

  // slice velocity data:
  //   --vel_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/dycoms_velocity.h5' \
  //   --init_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/dycoms_init.h5' \
  // piggy velocity data:
  //   --vel_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/velocity_out.dat' \

  // 1-mom micro opts:
  //   --cond=1 --cevp=1 --revp=1 --conv=1 --accr=1 --sedi=1 \
  // 2-mom micro opts:
  //   --acti=1 --cond=1 --accr=1 --acnv=1 --sedi=1 \

  //TODO - blk_1m piggy doesnt work with latent heat flux = 93 W/m2 (works ok for 50 W/m2) (in general it's in the rv_src)
  //TODO - small negative numbers in both schemes...

  string opts_common = 
    "--outfreq=200 --nt=1800 --spinup=0 --nx=97 --nz=301 --dt=1 --relax_th_rv=false";
  set<string> opts_micro({
    "--micro=blk_1m --outdir=out_blk_1m_piggy_7200 --adv_serial=false --async=true --backend=OpenMP --case=dycoms \
     --slice=false --piggy=true \
     --init_type='dat' \
     --init_dir='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/blk_1m_7200/'\
     --vel_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/blk_1m_7200/velocity_out.dat' \
     --w_src=0 --uv_src=0 --rv_src=1 --th_src=1 --subsidence=1 "
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
