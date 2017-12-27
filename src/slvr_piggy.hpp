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
      std::cout << "piggybacking from: " << vel_in << std::endl;

      slice = vm["slice"].as<bool>();
      std::cout << "slice flag: " << slice << std::endl;

      in_bfr.resize(this->state(this->vip_ixs[0]).shape());

      // open file for in vel
      // TODO: somehow check dimensionality of the input arrays
      if(!slice)
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
      else if(slice)
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
          f_vel_in >> in_bfr;
          this->state(this->vip_ixs[d]) = in_bfr;
        }
      }
      else if(slice)
      {
        std::cerr<<"slvr_piggy: read-in PyCLES velocity field"<<std::endl;
        using real_t = typename ct_params_t::real_t;
        using namespace libmpdataxx::arakawa_c;

        std::string fname  = vel_in + std::to_string(this->timestep)+".hdf";

        std::cerr<<"trying to read from file: "<<fname<<std::endl;

        H5::H5File h5f(fname, H5F_ACC_RDONLY);
        H5::Group h5g = h5f.openGroup("data_0");
        {
          // get horizontal velocity data set
          H5::Group h5gg  = h5g.openGroup("v");
          H5::DataSet h5d = h5gg.openDataSet("data_0");
          // ... and its dimension    
          H5::DataSpace h5s = h5d.getSpace();
          hsize_t data_dim[2];
          h5s.getSimpleExtentDims(data_dim, NULL);
          //data_dim[0]=4;
          //data_dim[1]=9;

          //create temporary arrays TODO - get rid of it
          blitz::Array<float, 2> tmp_les_v(data_dim[0],   data_dim[1]  );
          blitz::Array<float, 2>     tmp_v(data_dim[0]+1, data_dim[1]+1);

          // read horizontal velocity to a temporary array
          h5d.read(tmp_les_v.data(), H5::PredType::NATIVE_FLOAT);

          //tmp_les_v = 1.;
          //tmp_v     = 2.;

          blitz::Range row(1, data_dim[0]);
          blitz::Range col(1, data_dim[1]-1);
          blitz::Range col_all(0, data_dim[1]);
          blitz::Range row_all(0, data_dim[0]);

          tmp_v(row, col) = real_t(0.5) * (tmp_les_v(row-1, col-1) + tmp_les_v(row-1, col));
          tmp_v(row, 0) = tmp_les_v(row-1, 0);
          tmp_v(row, data_dim[1]) = tmp_les_v(row-1, data_dim[1]-1);
          tmp_v(0, col_all) = tmp_v(data_dim[0], col_all); //cyclic

          this->state(ix::vip_i) = 4.;
          for (int i=0; i<=data_dim[0]; i++){
              for(int j=0; j<=data_dim[1]; j++){
                  this->state(this->vip_ixs[0])(i,j) = tmp_v(i,j);
              }
          }
          //this->state(ix::vip_i)(row_all, col_all) = tmp_v(row_all, col_all); TODO - why is this not working?

          // fill halo
          this->xchng_sclr(this->state(this->vip_ixs[0]), this->ijk, this->halo);

          tmp_les_v.free();
          tmp_v.free();
        }
        {
          // get vertical velocity data set
          H5::Group h5gg  = h5g.openGroup("w");
          H5::DataSet h5d = h5gg.openDataSet("data_0");
          // ... and its dimension    
          H5::DataSpace h5s = h5d.getSpace();
          hsize_t data_dim[2];
          h5s.getSimpleExtentDims(data_dim, NULL);
          //data_dim[0]=4;
          //data_dim[1]=9;

          //create temporary arrays TODO - get rid of it
          blitz::Array<float, 2> tmp_les_w(data_dim[0],   data_dim[1]  );
          blitz::Array<float, 2>     tmp_w(data_dim[0]+1, data_dim[1]+1);

          // read horizontal velocity to a temporary array
          h5d.read(tmp_les_w.data(), H5::PredType::NATIVE_FLOAT);

          //tmp_les_w = 1.;
          //tmp_w     = 2.;

          blitz::Range row(1, data_dim[0]-1);
          blitz::Range col(1, data_dim[1]);
          blitz::Range col_all(0, data_dim[1]);
          blitz::Range row_all(0, data_dim[0]);

          tmp_w(row, col) = real_t(0.5) * (tmp_les_w(row-1, col-1) + tmp_les_w(row, col-1));
          tmp_w(0, col)   = real_t(0.5) * (tmp_les_w(0, col-1) + tmp_les_w(data_dim[0]-1, col-1));
          tmp_w(data_dim[0], col) = tmp_w(0, col); //cyclic
          tmp_w(row_all, 0) = 0.;                  //wall

          for (int i=0; i<=data_dim[0]; i++){
              for(int j=0; j<=data_dim[1]; j++){
                  this->state(this->vip_ixs[1])(i,j) = tmp_w(i,j);
              }
          }
          //this->state(ix::vip_j)(row_all, col_all) = tmp_w(row_all, col_all); TODO - why is this not working?

          // fill halo
          this->xchng_sclr(this->state(this->vip_ixs[1]), this->ijk, this->halo);

          tmp_les_w.free();
          tmp_w.free();
        }
      }
      else
      {
        assert(false);
      }
    }
    //std::cerr<<"vip_ixs[0] = " << this->state(this->vip_ixs[0]) <<std::endl;
    //std::cerr<<"vip_ixs[1] = " << this->state(this->vip_ixs[1]) <<std::endl;

    this->mem->barrier();
  }

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p) {}
};

