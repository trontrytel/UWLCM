#pragma once
#include <libmpdata++/solvers/mpdata_rhs.hpp>
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
  typename std::enable_if<ct_params_t::piggy == 1 && ct_params_t::slice == 0>::type
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
      //in_bfr.resize(this->state(this->vip_ixs[0]).shape());

      // open file for in vel
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
 
      for (int d = 0; d < parent_t::n_dims; ++d)
      {

        in_bfr.resize(this->state(this->vip_ixs[d]).shape());
        // read in through buffer, if done directly caused data races
        // TODO - change to hdf5?
        f_vel_in >> in_bfr;
        this->state(this->vip_ixs[d]) = in_bfr;
        in_bfr.resize(0);
if (d==0)
{
std::cerr<<"reading " << vel_in << " timestep "<< std::to_string(this->timestep)<<std::endl;
std::cerr<<" "<<std::endl;
std::cerr<<"rv (min, max) = (" << blitz::min(this->state(ix::rv)) << " , " << blitz::max(this->state(ix::rv)) << ")" << std::endl;
std::cerr<<"th (min, max) = (" << blitz::min(this->state(ix::th)) << " , " << blitz::max(this->state(ix::th)) << ")" << std::endl;
}
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

// slice
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 1 && ct_params_t::slice == 1>::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs<ct_params_t>
  >
{
  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs<ct_params_t>
  >;  

  private:
  typename parent_t::arr_t tmp_v, tmp_w, tmp_div; // input buffer for velocity
  H5::Group h5g_v, h5g_w;                         // hdf5 group, dataset and dimensions of velocity data
  H5::DataSet h5d_v, h5d_w;              
  H5::DataSpace h5s_v, h5s_w;
  hsize_t data_dim_v[2], data_dim_w[2];
 
  protected:
  std::ifstream f_vel_in; // input velocity file
  std::string vel_in;

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      po::options_description opts("Slice options"); 
      opts.add_options()
        ("vel_in", po::value<std::string>()->required(), "file with input velocities")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
      vel_in = vm["vel_in"].as<std::string>();

      // open file for in vel
      // TODO check dimensionality of the input arrays
      try
      {
        std::cout << "piggybacking from: " << vel_in << std::endl;
        H5::H5File h5f(vel_in, H5F_ACC_RDONLY);

        // get velocity data sets ...
        h5g_v   = h5f.openGroup("v_nodiv_dual"); // TODO - move this to ante loop
        h5g_w   = h5f.openGroup("w_nodiv_dual");

        // ... the correct time step (t=0) ...
        h5d_v = h5g_v.openDataSet(std::to_string(0));
        h5d_w = h5g_w.openDataSet(std::to_string(0));

        // ... and data dimensions
        h5s_v = h5d_v.getSpace();
        h5s_w = h5d_w.getSpace();
        h5s_v.getSimpleExtentDims(data_dim_v, NULL);
        h5s_w.getSimpleExtentDims(data_dim_w, NULL);
        tmp_v.resize(data_dim_v[0], data_dim_v[1]);
        tmp_w.resize(data_dim_w[0], data_dim_w[1]);
        tmp_div.resize(data_dim_w[0], data_dim_v[1]);
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

    // read velocity in advector
    //TODO - should be multiplied by rho
    if(this->rank==0)
    {
      using ix = typename ct_params_t::ix;
      using namespace libmpdataxx::arakawa_c;
      std::cerr<<"reading " << vel_in << " timestep "<< std::to_string(this->timestep)<<std::endl;

      // ... the correct time step ...
      h5d_v = h5g_v.openDataSet(std::to_string(this->timestep));
      h5d_w = h5g_w.openDataSet(std::to_string(this->timestep));

      h5d_v.read(tmp_v.data(), H5::PredType::NATIVE_FLOAT);
      h5d_w.read(tmp_w.data(), H5::PredType::NATIVE_FLOAT);

      // read the data from temporary array to vip array
      for (int i=0; i<data_dim_v[0]; i++){
          for(int j=0; j<data_dim_v[1]; j++){
              this->mem->GC[0](i,j) = tmp_v(i,j) / 35.; //* //(this->mem->g_factor()(i,j) + this->mem->g_factor()(i+1,j)) * 2;   //TODO
          }
      }
      // read the data from temporary array to vip array
      for (int i=0; i<data_dim_w[0]; i++){
          for(int j=0; j<data_dim_w[1]; j++){
              this->mem->GC[1](i,j) = tmp_w(i,j) / 5.;// * //(this->mem->g_factor()(i,j) + this->mem->g_factor()(i, j+1)) * 2;   //TODO
          }
      }

      // fill halo
      this->xchng_vctr_alng(this->mem->GC, true);
      this->xchng_vctr_nrml(this->mem->GC, this->ijk);

      //typename ct_params_t::real_t max_abs_div_eps = 1e-3;
      for (int i=0; i<data_dim_w[0]; i++){
          for(int j=0; j<data_dim_w[1]; j++){
            tmp_div(i,j) = (
                      (this->mem->GC[0](i+1, j) - this->mem->GC[0](i, j)) + 
                      (this->mem->GC[1](i, j+1) - this->mem->GC[1](i, j))
                     ) / this->mem->g_factor()(i, j);
          }
      }
      typename ct_params_t::real_t max_abs_div = max(abs(tmp_div));
      std::cerr<<"max abs div = "<<max_abs_div<<std::endl;

      //if (max_abs_div > this->max_abs_div_eps)
      //    throw std::runtime_error("initial advector field is divergent");

std::cerr<<" "<<std::endl;
std::cerr<<"rv (min, max) = (" << blitz::min(this->state(ix::rv)) << " , " << blitz::max(this->state(ix::rv)) << ")" << std::endl;
std::cerr<<"th (min, max) = (" << blitz::min(this->state(ix::th)) << " , " << blitz::max(this->state(ix::th)) << ")" << std::endl;
std::cerr<<"-------------------------------------------"<<std::endl;
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
