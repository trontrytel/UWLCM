#pragma once
#include "slvr_blk_1m_common.hpp"

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t, class enableif = void>
class slvr_blk_1m 
{};

// 2D version 
template <class ct_params_t>
class slvr_blk_1m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_1m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_1m( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) : 
    parent_t(args, p)
  {}  
  
  protected:
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

// 3D version 
template <class ct_params_t>
class slvr_blk_1m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_blk_1m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_1m( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) : 
    parent_t(args, p)
  {}  

  protected:
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        auto 
  	dot_rr = rhs.at(parent_t::ix::rr)(i, j, this->k);
        const auto 
        rhod   = (*this->mem->G)(i, j, this->k),
  	rr     = this->state(parent_t::ix::rr)(i, j, this->k);
        libcloudphxx::blk_1m::rhs_columnwise<real_t>(this->opts, dot_rr, rhod, rr, this->params.dz);
      }
  }
};
