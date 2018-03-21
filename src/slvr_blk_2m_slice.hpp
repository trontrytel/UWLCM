#pragma once
#include "slvr_blk_2m.hpp"

template <class ct_params_t, class enableif = void>
class slvr_blk_2m_slice 
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx;

// 2D version 
template <class ct_params_t>
class slvr_blk_2m_slice<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_2m<ct_params_t>
{
  public:
  using parent_t = slvr_blk_2m<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_2m_slice( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) : 
    parent_t(args, p)
  {} 

  protected:


  void hook_post_step()
  {
    using ix = typename ct_params_t::ix;
    using namespace libmpdataxx::arakawa_c;
    for(auto a: std::list<int>({ix::rc, ix::rr, ix::rv, ix::th, ix::nc, ix::nr}))
    {
      //this->state(a)(this->ijk) /= this->state(ix::one)(this->ijk);
      this->xchng(a);
    }
 
    //this->state(ix::one)(this->ijk) /= this->state(ix::one)(this->ijk);
    //this->xchng(ix::one);
 
    parent_t::hook_post_step(); // includes forcings

    //TODO - some problem with surface fluxes
    this->state(ix::rc)(this->i, 0) = 0.;

    //std::cerr<<"rc (min, max) = (" << blitz::min(this->state(ix::rc)) << " , " << blitz::max(this->state(ix::rc)) << ")" << std::endl;
    //std::cerr<<"nc (min, max) = (" << blitz::min(this->state(ix::nc)) << " , " << blitz::max(this->state(ix::nc)) << ")" << std::endl;

    //TODO - TMP!!!
    this->state(ix::nc)(this->ijk) = blitz::where(
      this->state(ix::nc)(this->ijk) < 0.,
      0.,
      this->state(ix::nc)(this->ijk)
    );
  }
};
