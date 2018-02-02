#pragma once
#include "slvr_blk_2m_common.hpp"

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t, class enableif = void>
class slvr_blk_2m 
{};

// 2D version 
template <class ct_params_t>
class slvr_blk_2m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_2m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_2m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_2m( 
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

    using ix = typename parent_t::ix;
    // column-wise

    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      auto
        dot_rr = rhs.at(ix::rr)(i, this->j),
        dot_nr = rhs.at(ix::nr)(i, this->j);
      const auto
        rhod   = (*this->mem->G)(i, this->j),
        rr     = this->state(ix::rr)(i, this->j),
        nr     = this->state(ix::nr)(i, this->j);

      libcloudphxx::blk_2m::rhs_columnwise<real_t>(
        this->opts, dot_rr, dot_nr,
        rhod,     rr,     nr,
        this->dt,
        this->params.dz
      );
    }
  }
};

