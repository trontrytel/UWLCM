#pragma once
#include "slvr_blk_1m_common.hpp"

template <class ct_params_t, class enableif = void>
class slvr_blk_1m_slice 
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

// 2D version 
template <class ct_params_t>
class slvr_blk_1m_slice<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_1m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_1m_slice( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) : 
    parent_t(args, p)
  {} 

  protected:

  void hook_ante_loop(int nt)
  {
    using ix = typename ct_params_t::ix;
    this->state(ix::one)(this->i, this->j) = 1.;

    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }


  void hook_post_step()
  {
    //TODO
    using ix = typename ct_params_t::ix;
    using namespace libmpdataxx::arakawa_c;
    for(auto a: std::list<int>({ix::rc, ix::rr, ix::rv, ix::th}))
    {
      this->state(a)(this->i, this->j) /= this->state(ix::one)(this->i, this->j);
      //this->xchng_sclr(this->state(a), this->i^this->halo, this->j^this->halo);
      this->xchng(a);
    }
    this->state(ix::one)(this->i, this->j) /= this->state(ix::one)(this->i, this->j);
    this->xchng(ix::one);
 
    parent_t::condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings

   this->state(ix::rc)(this->i, 0) = 0.; //TODO
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?
    
    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    { 
      auto 
        dot_rr = rhs.at(parent_t::ix::rr)(i, this->j);
      const auto 
        rhod   = (*this->mem->G)(i, this->j),
        rr     = this->state(parent_t::ix::rr)(i, this->j);
      libcloudphxx::blk_1m::rhs_columnwise<real_t>(this->opts, dot_rr, rhod, rr, this->params.dz);
    }
  }
};
