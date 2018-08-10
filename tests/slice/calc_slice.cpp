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

  //  "--outfreq=200 --nt=12000 --spinup=9600 --nx=97 --nz=301 --dt=0.75 --relax_th_rv=false";
  //  "--outfreq=200 --nt=2400 --spinup=0 --nx=97 --nz=301 --dt=0.75 --relax_th_rv=false";
  //   --init_dir='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/blk_2m_9600/'\
  //   --vel_in='/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/blk_2m_9600/velocity_out.dat' \

  //--acnv_A=1350 --acnv_b=2.47 --acnv_c=-1.79 \
  //--blk2m_mean_rd=0.05e-6 --blk2m_sdev_rd=1.5 --blk2m_N_stp=080e6 \
  //  --acnv_A=1350  --acnv_b=2.47  --acnv_c=-1.79  2-mom
  //  --rc_0=5e-4 --k_acnv=1e-3"                    1-mom

  string opts_common = 
    "--outfreq=200 --nt=9000 --spinup=7200 --nx=97 --nz=301 --dt=1 --relax_th_rv=false";
  set<string> opts_micro({
    "--micro=lgrngn --outdir=out_lgrngn --adv_serial=false --async=true --backend=CUDA --case=dycoms \
     --slice=false --piggy=false \
     --init_type=calc \
     --init_dir='none'\
     --vel_in='none' \
     --sd_conc=32 --sstp_cond=1 --z_rlx_sclr=100 --sstp_coal=1\
     --w_src=1 --uv_src=1 --rv_src=1 --th_src=1 --subsidence=1 "
      " --out_wet=\""
        ".5e-6:25e-6|0,1,2,3,6;"  // FSSP
        "25e-6:1|0,1,2,3,6;"      // "rain"
        ".5e-6:1|0,1,2,3,6;"      // all hydro
        + bins_wet_str + // aerosol spectrum (wet)
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
