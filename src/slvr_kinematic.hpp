#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include "detail/checknan.cpp"

using namespace libmpdataxx; // TODO: get rid of it?

// piggybacker
template <class ct_params_t>
class slvr_kinematic<
  ct_params_t
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip<ct_params_t>
  >
{
  protected:

  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip<ct_params_t>
  >;  

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      po::options_description opts("Kinematic options"); 
      opts.add_options()
        ("vel_in", po::value<std::string>()->required(), "file with input velocities")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
          
      vel_in = vm["vel_in"].as<std::string>();
      std::cout << "piggybacking from: " << vel_in << std::endl;

      in_bfr.resize(this->state(this->vip_ixs[0]).shape());
      // open file for in vel
      // TODO: somehow check dimensionality of the input arrays
      try{
        f_vel_in.open(vel_in); 
      }
      catch(...)
      {
        throw std::runtime_error("error opening velocities input file defined by --vel_in");
      }
      this->record_aux_const("piggybacking", -44); // dummy -44 
      this->record_aux_const(std::string("vel_in : ") + vel_in, -44);  // dummy -44
    }
    this->mem->barrier();
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // do whatever
    this->mem->barrier(); //necessary?
    // read velo, overwrite any vel rhs
    if(this->rank==0)
    {
      using ix = typename ct_params_t::ix;

      for (int d = 0; d < parent_t::n_dims; ++d)
      {
        // read in through buffer, if done directly caused data races
        f_vel_in >> in_bfr;
        this->state(this->vip_ixs[d]) = in_bfr;
//std::cout << this->state(this->vip_ixs[d]);
      }
    }
    this->mem->barrier();
  }

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p) {}
};

