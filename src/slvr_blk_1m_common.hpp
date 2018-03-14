#pragma once
#include "slvr_common.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

template <class ct_params_t>
class slvr_blk_1m_common : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix; // TODO: it's now in solver_common - is it needed here?
  using real_t = typename ct_params_t::real_t;
  protected:

  void condevap()
  {
    auto 
      th   = this->state(ix::th)(this->ijk), // potential temperature
      rv   = this->state(ix::rv)(this->ijk), // water vapour mixing ratio
      rc   = this->state(ix::rc)(this->ijk), // cloud water mixing ratio
      rr   = this->state(ix::rr)(this->ijk); // rain water mixing ratio
    auto const
      rhod = (*this->mem->G)(this->ijk);

   
    libcloudphxx::blk_1m::adj_cellwise<real_t>( 
      opts, rhod, th, rv, rc, rr, this->dt
    );
    this->mem->barrier(); 
  }

  void zero_if_uninitialised(int e)
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
      this->state(e)(this->ijk) = 0;
  }

  bool get_rain() { return opts.conv; }
  void set_rain(bool val) { opts.conv = val; };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rc);
    zero_if_uninitialised(ix::rr);

    if (this->params.init_type == "dat")
    {
      std::ifstream rc_in;
      blitz::Array<double, 2> in_bfr; //has to be double to properly read in nc and rc data
     
      in_bfr.resize(this->state(ix::rc).shape());
      rc_in.open(this->params.init_dir+"rc.dat");
      rc_in >> in_bfr;
//TODO - it should be enough to change 2 to something better
//TODO - when fixed, uncomment the 3D version in bicycles
///Users/ajaruga/clones/UWLCM/src/slvr_blk_1m_common.hpp:62:27: 
//note: in instantiation of function template specialization 
//'blitz::Array<float, 3>::operator = <blitz::Array<double, 2> >' requested here
//this->state(ix::rc) = in_bfr;

      this->state(ix::rc) = in_bfr;

      rc_in.close();
    }

    // deal with initial supersaturation
    condevap();

    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

  void hook_ante_step()
  {
    parent_t::hook_ante_step();
    // store rl for buoyancy
    this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    // cell-wise
    {
      auto 
	dot_rc = rhs.at(ix::rc)(this->ijk),
	dot_rr = rhs.at(ix::rr)(this->ijk);
      const auto 
	rc   = this->state(ix::rc)(this->ijk),
	rr   = this->state(ix::rr)(this->ijk);
      libcloudphxx::blk_1m::rhs_cellwise<real_t>(opts, dot_rc, dot_rr, rc, rr);
    }
  }

  // 
  void hook_post_step()
  {
    //TODO TMP should be removed after I read in divergence free velocity data
    this->cleanup(ix::rv);
    this->cleanup(ix::rc);
    this->cleanup(ix::rr);

if (this->rank == 0)
{
std::cerr<<" "<<std::endl;
std::cerr<<"timestep = "<<this->timestep<<std::endl;
std::cerr<<"rv (min, max) = (" << blitz::min(this->state(ix::rv)) << " , " << blitz::max(this->state(ix::rv)) << ")" << std::endl;
std::cerr<<"rc (min, max) = (" << blitz::min(this->state(ix::rc)) << " , " << blitz::max(this->state(ix::rc)) << ")" << std::endl;
std::cerr<<"rr (min, max) = (" << blitz::min(this->state(ix::rr)) << " , " << blitz::max(this->state(ix::rr)) << ")" << std::endl;
std::cerr<<" "<<std::endl;
}

    condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment

/*
if (this->timestep == 7200 && this->rank == 0){
std::ofstream th_out_init, rv_out_init, rc_out_init, rr_out_init;
th_out_init.open(this->outdir+"/th_out_init_7200.dat");
rv_out_init.open(this->outdir+"/rv_out_init_7200.dat");
rc_out_init.open(this->outdir+"/rc_out_init_7200.dat");
rr_out_init.open(this->outdir+"/rr_out_init_7200.dat");
th_out_init << this->state(ix::th)(this->ijk);
rv_out_init << this->state(ix::rv)(this->ijk);
rc_out_init << this->state(ix::rc);
rr_out_init << this->state(ix::rr)(this->ijk);
th_out_init.close();
rv_out_init.close();
rc_out_init.close();
rr_out_init.close();
}
*/

    parent_t::hook_post_step(); // includes the above forcings

    this->mem->barrier();
  }

  libcloudphxx::blk_1m::opts_t<real_t> opts;

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::blk_1m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  slvr_blk_1m_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  {}  
};

