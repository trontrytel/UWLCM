#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include "detail/checknan.cpp"

template <class ct_params_t, class enableif = void>
class slvr_piggy
{};

using namespace libmpdataxx; // TODO: get rid of it?

// driver
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 0 >::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs<ct_params_t>
  >
{
  private:
  bool save_vel; // should velocity field be stored for piggybacking

  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs<ct_params_t>
  >;  

  std::ofstream f_vel_out; // file for velocity field

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      po::options_description opts("Driver options"); 
      opts.add_options()
        ("save_vel", po::value<bool>()->default_value(false), "should velocity field be stored for piggybacking")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
          
      save_vel = vm["save_vel"].as<bool>();
      // open file for out vel
      if(save_vel)
      {
        try{
          f_vel_out.open(this->outdir+"/velocity_out.dat"); 
        }
        catch(...)
        {
          throw std::runtime_error("error opening velocity output file '{outdir}/velocity_out.dat'");
        }
      }
    }
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes changes of velocity field due to vip_rhs_impl_fnlz()
    this->mem->barrier();
    // save velocity field
    if(this->rank==0 && save_vel)
    {
      for (int d = 0; d < parent_t::n_dims; ++d)
      {
        f_vel_out << this->state(this->vip_ixs[d]);
      }
    }
  }

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t p
  ) :
    parent_t(args, p) {}

  public:

  // ctor
  struct rt_params_t : parent_t::rt_params_t 
  {   
    rt_params_t()
    {
      this->prs_tol = 1e-6;
    }
  }; 
};


// piggybacker
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 1 >::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip<ct_params_t>
  >
{

  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip<ct_params_t>
  >;  

  private:
  typename parent_t::arr_t in_bfr; // input buffer for velocity
  std::string vel_in;
  
  protected:

  std::ifstream f_vel_in; // input velocity file

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      po::options_description opts("Piggybacker options"); 
      opts.add_options()
        ("vel_in", po::value<std::string>()->required(), "file with input velocities")
        ("pycles", po::value<bool>()->default_value(false), "is velocity field from PyCLES")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
          
      vel_in = vm["vel_in"].as<std::string>();
      std::cout << "piggybacking from: " << vel_in << std::endl;

      user_params_t user_params;
      user_params.pycles = vm["pycles"].as<bool>();

      in_bfr.resize(this->state(this->vip_ixs[0]).shape());

      // open file for in vel
      // TODO: somehow check dimensionality of the input arrays
      if(!user_params.pycles)
      {
        try
        {
          f_vel_in.open(vel_in);
        }
        catch(...)
        {
          throw std::runtime_error("error opening velocities input file defined by --vel_in");
        }
      }
      else if(user_params.pycles)
      {  
        try
        {
          f_vel_in.open(vel_in + "1.hdf"); 
        }
        catch(...)
        {
          throw std::runtime_error("error opening velocities input file defined by --vel_in");
        }
      }
      else
      {
        assert(false);
      }
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

      //TODO - how to pass it from user-defined options?
      bool pycles = true;
 
      if(!pycles)
      {
        for (int d = 0; d < parent_t::n_dims; ++d)
        {
          // read in through buffer, if done directly caused data races
          f_vel_in >> in_bfr;
          this->state(this->vip_ixs[d]) = in_bfr;
          //std::cout << this->state(this->vip_ixs[d]);
        }
      }
      else if(pycles)
      {
        for (int d = 0; d < parent_t::n_dims; ++d)
        {
          std::cerr<<"AQQ - here we should read my velocity"<<std::endl;
          this->state(this->vip_ixs[d]) = 0.;
        }
        this->mem->advectee(ix::u) = 0.;
        this->mem->advectee(ix::w) = 0.;
        this->mem->GC[0] = 0;
        this->mem->GC[1] = 0;
        this->mem->advector(ix::u) = 0.;
        this->mem->advector(ix::w) = 0.;
      }
      else
      {
        assert(false);
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

