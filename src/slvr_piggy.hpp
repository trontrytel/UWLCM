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
  
  protected:

  std::ifstream f_vel_in; // input velocity file
  std::string vel_in;
  bool slice;

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      po::options_description opts("Piggybacker options"); 
      opts.add_options()
        ("vel_in", po::value<std::string>()->required(), "file with input velocities")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
          
      vel_in = vm["vel_in"].as<std::string>();
      slice  = vm["slice"].as<bool>();

      std::cout << "piggybacking from: " << vel_in << std::endl;

      in_bfr.resize(this->state(this->vip_ixs[0]).shape());

      // open file for in vel
      // TODO: somehow check dimensionality of the input arrays
      // TODO: check if this works for hdf5 files
      try
      {
        f_vel_in.open(vel_in);
      }
      catch(...)
      {
        throw std::runtime_error("error opening velocities input file defined by --vel_in");
      }
    }
    this->mem->barrier();
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // do whatever
    this->mem->barrier();

    // read velo, overwrite any vel rhs
    if(this->rank==0)
    {
      using ix = typename ct_params_t::ix;

      if(!slice)
      {
        for (int d = 0; d < parent_t::n_dims; ++d)
        {
          // read in through buffer, if done directly caused data races
          // TODO - change to hdf5 
          f_vel_in >> in_bfr;
          this->state(this->vip_ixs[d]) = in_bfr;
        }
      }
      else if(slice)
      {
        using namespace libmpdataxx::arakawa_c;

        std::string fname  = vel_in; 
        std::cerr<<"reading "<<fname<<" timestep "<< std::to_string(this->timestep)<<std::endl;
        H5::H5File h5f(fname, H5F_ACC_RDONLY);

        // get velocity data sets ...
        //H5::Group h5g_v   = h5f.openGroup("v"); // TODO - move this to ante loop
        //H5::Group h5g_w   = h5f.openGroup("w");
        H5::Group h5g_v   = h5f.openGroup("v_nodiv"); // TODO - move this to ante loop
        H5::Group h5g_w   = h5f.openGroup("w_nodiv");
        // ... the correct time step ...
        H5::DataSet h5d_v = h5g_v.openDataSet(std::to_string(this->timestep));
        H5::DataSet h5d_w = h5g_w.openDataSet(std::to_string(this->timestep));
        // ... and data dimensions
        H5::DataSpace h5s = h5d_v.getSpace();
        hsize_t data_dim[2];
        h5s.getSimpleExtentDims(data_dim, NULL);

        // create temporary array and read the data 
        // TODO - can I read in directly?
        // TODO - use the same bufer as for non slice? are the dimensions the same for bufer or are they bigger because they include halo?
        blitz::Array<float, 2> tmp_v(data_dim[0], data_dim[1]);
        blitz::Array<float, 2> tmp_w(data_dim[0], data_dim[1]);
        h5d_v.read(tmp_v.data(), H5::PredType::NATIVE_FLOAT);
        h5d_w.read(tmp_w.data(), H5::PredType::NATIVE_FLOAT);

std::cerr<<"-------------------------------------------"<<std::endl;
std::cerr<<"tmp_v (min, max) = (" << blitz::min(tmp_v) << " , " << blitz::max(tmp_v) << ")" << std::endl;
std::cerr<<"tmp_w (min, max) = (" << blitz::min(tmp_w) << " , " << blitz::max(tmp_w) << ")" << std::endl;

        //this->state(this->vip_ixs[0])(this->ijk) = 0.;
        //this->state(this->vip_ixs[1])(this->ijk) = 0.;
 
        // read the data from temporary array to vip array
        for (int i=0; i<data_dim[0]; i++){
            for(int j=0; j<data_dim[1]; j++){
                this->state(this->vip_ixs[0])(i,j) = tmp_v(i,j);
                this->state(this->vip_ixs[1])(i,j) = tmp_w(i,j);
            }
        }

        // fill halo
        this->xchng_sclr(this->state(this->vip_ixs[0]), this->ijk, this->halo);
        this->xchng_sclr(this->state(this->vip_ixs[1]), this->ijk, this->halo);

std::cerr<<"uwlcm_v (min, max) = (" << blitz::min(this->state(this->vip_ixs[0])) << " , " << blitz::max(this->state(this->vip_ixs[0])) << ")" << std::endl;
std::cerr<<"uwlcm_w (min, max) = (" << blitz::min(this->state(this->vip_ixs[1])) << " , " << blitz::max(this->state(this->vip_ixs[1])) << ")" << std::endl;

std::cerr<<" "<<std::endl;
std::cerr<<"th (min, max) = (" << blitz::min(this->state(ix::th)) << " , " << blitz::max(this->state(ix::th)) << ")" << std::endl;
std::cerr<<"rv (min, max) = (" << blitz::min(this->state(ix::rv)) << " , " << blitz::max(this->state(ix::rv)) << ")" << std::endl;
//std::cerr<<"rc (min, max) = (" << blitz::min(this->state(ix::rc)) << " , " << blitz::max(this->state(ix::rc)) << ")" << std::endl;
//std::cerr<<"nc (min, max) = (" << blitz::min(this->state(ix::nc)) << " , " << blitz::max(this->state(ix::nc)) << ")" << std::endl;

std::cerr<<"-------------------------------------------"<<std::endl;

        // cleanup
        tmp_v.free();
        tmp_w.free();
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

