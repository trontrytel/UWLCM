#pragma once
#include "slvr_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

template <class ct_params_t>
class slvr_blk_2m_common : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix; // TODO: it's now in solver_common - is it needed here?
  using real_t = typename ct_params_t::real_t;
  protected:

  void zero_if_uninitialised(int e)
  {
    if (!finite(sum(this->state(e)(this->ijk))))
      this->state(e)(this->ijk) = 0;
  }

  bool get_rain() { return opts.acnv; }
  void set_rain(bool val) 
  { 
    opts.acnv = val; 
    opts.RH_max = val ? 44 : 1.01;
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    if (this->params.init_type == "dat")
    { 
      std::ifstream nc_in, rc_in;
      blitz::Array<double, 2> in_bfr; //has to be double to properly read in nc and rc data
      in_bfr.resize(this->state(ix::nc).shape());

      rc_in.open(this->params.init_dir+"rc.dat");
      nc_in.open(this->params.init_dir+"nc.dat");

      nc_in >> in_bfr;
      this->state(ix::nc) = in_bfr;
      rc_in >> in_bfr;
      this->state(ix::rc) = in_bfr;

      rc_in.close();
      nc_in.close();
    }

    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rc);
    zero_if_uninitialised(ix::nc);
    zero_if_uninitialised(ix::rr);
    zero_if_uninitialised(ix::nr);

    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

  void hook_ante_step()
  {
    this->mem->barrier();
    parent_t::hook_ante_step();

    // store rl for buoyancy
    this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);
  }

  void hook_post_step()
  {
    //TODO TMP - should be removed after I read in divergence free files 
    //from my Helmholtz decomposition 
    this->cleanup(ix::nc); //<- the most important one
    this->cleanup(ix::nr);
    this->cleanup(ix::rc);
    this->cleanup(ix::rr);
    this->cleanup(ix::rv);
    /*
if (this->rank==0){
std::cerr<<" "<<std::endl;
std::cerr<<"timestep = "<<this->timestep<<std::endl;
std::cerr<<"nc (min, max) = (" << blitz::min(this->state(ix::nc)) << " , " << blitz::max(this->state(ix::nc)) << ")" << std::endl;
std::cerr<<"rc (min, max) = (" << blitz::min(this->state(ix::rc)) << " , " << blitz::max(this->state(ix::rc)) << ")" << std::endl;
std::cerr<<"nr (min, max) = (" << blitz::min(this->state(ix::nr)) << " , " << blitz::max(this->state(ix::nr)) << ")" << std::endl;
std::cerr<<"rr (min, max) = (" << blitz::min(this->state(ix::rr)) << " , " << blitz::max(this->state(ix::rr)) << ")" << std::endl;
std::cerr<<"rv (min, max) = (" << blitz::min(this->state(ix::rv)) << " , " << blitz::max(this->state(ix::rv)) << ")" << std::endl;
std::cerr<<"th (min, max) = (" << blitz::min(this->state(ix::th)) << " , " << blitz::max(this->state(ix::th)) << ")" << std::endl;
std::cerr<<" "<<std::endl;
}
    */
/*
if (this->timestep == 7200 && this->rank == 0){
std::ofstream th_out_init, rv_out_init, nc_out_init, rc_out_init, nr_out_init, rr_out_init;
th_out_init.open(this->outdir+"/th_out_init_7200.dat");
rv_out_init.open(this->outdir+"/rv_out_init_7200.dat");
nc_out_init.open(this->outdir+"/nc_out_init_7200.dat");
rc_out_init.open(this->outdir+"/rc_out_init_7200.dat");
nr_out_init.open(this->outdir+"/nr_out_init_7200.dat");
rr_out_init.open(this->outdir+"/rr_out_init_7200.dat");
th_out_init << this->state(ix::th)(this->ijk);
rv_out_init << this->state(ix::rv)(this->ijk);
nc_out_init << this->state(ix::nc)(this->ijk);
rc_out_init << this->state(ix::rc)(this->ijk);
nr_out_init << this->state(ix::nr)(this->ijk);
rr_out_init << this->state(ix::rr)(this->ijk);
th_out_init.close();
rv_out_init.close();
nc_out_init.close();
rc_out_init.close();
nr_out_init.close();
rr_out_init.close();
}
*/
    parent_t::hook_post_step();
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at);

    this->mem->barrier(); // TODO: if neccesarry, then move to adv_rhs/....hpp

    // cell-wise
    {
      auto
        dot_th = rhs.at(ix::th)(this->ijk),
        dot_rv = rhs.at(ix::rv)(this->ijk),
        dot_rc = rhs.at(ix::rc)(this->ijk),
        dot_rr = rhs.at(ix::rr)(this->ijk),
        dot_nc = rhs.at(ix::nc)(this->ijk),
        dot_nr = rhs.at(ix::nr)(this->ijk),
        rc     = this->state(ix::rc)(this->ijk),
        rr     = this->state(ix::rr)(this->ijk),
        nc     = this->state(ix::nc)(this->ijk),
        nr     = this->state(ix::nr)(this->ijk);

      const auto
        rhod   = (*this->mem->G)(this->ijk),
        th     = this->state(ix::th)(this->ijk),
        rv     = this->state(ix::rv)(this->ijk);

      libcloudphxx::blk_2m::rhs_cellwise<real_t>(
        opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
        rhod,     th,     rv,     rc,     nc,     rr,     nr,
        this->dt
      );
    }
    this->mem->barrier(); // TODO: if needed, move to adv+rhs
  }

  libcloudphxx::blk_2m::opts_t<real_t> opts;

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::blk_2m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  slvr_blk_2m_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  {}  
};

