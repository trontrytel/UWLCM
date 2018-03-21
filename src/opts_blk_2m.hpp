/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "opts_common.hpp"
#include "slvr_blk_2m_slice.hpp"

// simulation and output parameters for micro=blk_2m
template <class solver_t, class user_params_t, class case_ptr_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params, 
  const user_params_t &user_params,
  const case_ptr_t &case_ptr,
  typename std::enable_if<std::is_same<
    decltype(solver_t::rt_params_t::cloudph_opts),
    libcloudphxx::blk_2m::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Double-moment bulk microphysics options"); 
  opts.add_options()
    ("acti", po::value<bool>()->default_value(rt_params.cloudph_opts.acti) , "TODO (on/off)")
    ("cond", po::value<bool>()->default_value(rt_params.cloudph_opts.cond) , "TODO (on/off)")
    ("accr", po::value<bool>()->default_value(rt_params.cloudph_opts.accr) , "TODO (on/off)")
    ("acnv", po::value<bool>()->default_value(rt_params.cloudph_opts.acnv) , "TODO (on/off)")
    ("sedi", po::value<bool>()->default_value(rt_params.cloudph_opts.sedi) , "TODO (on/off)")
    ("acnv_A", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_A), "parameter in autoconversion rate formulae")
    ("acnv_b", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_b), "parameter in autoconversion rate formulae")
    ("acnv_c", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_c), "parameter in autoconversion rate formulae")
  ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Morrison and Grabowski 2007 scheme options
  rt_params.cloudph_opts.acti = vm["acti"].as<bool>();
  rt_params.cloudph_opts.cond = vm["cond"].as<bool>();
  rt_params.cloudph_opts.accr = vm["accr"].as<bool>();
  rt_params.cloudph_opts.acnv = vm["acnv"].as<bool>();
  rt_params.cloudph_opts.sedi = vm["sedi"].as<bool>();
  rt_params.cloudph_opts.acnv_A = vm["acnv_A"].as<typename solver_t::real_t>();
  rt_params.cloudph_opts.acnv_b = vm["acnv_b"].as<typename solver_t::real_t>();
  rt_params.cloudph_opts.acnv_c = vm["acnv_c"].as<typename solver_t::real_t>();

  rt_params.cloudph_opts.dry_distros.push_back({
    .mean_rd = 0.02e-6,
    .sdev_rd = 1.4,
    .N_stp   = 60e6,
    .chem_b  = .55
  });
/*
  rt_params.cloudph_opts.dry_distros.push_back({
    .mean_rd = setup.mean_rd2 / si::metres,
    .sdev_rd = setup.sdev_rd2,
    .N_stp   = setup.n2_stp * si::cubic_metres,
    .chem_b  = setup.chem_b
  });
*/

  // output variables
  rt_params.outvars = {
    // <TODO>: make it common among all three micro?
    {solver_t::ix::th, {"th", "[K]"}},
    {solver_t::ix::rv, {"rv", "[kg kg-1]"}},
    {solver_t::ix::w,  {"w", "[m/s]"}},
    {solver_t::ix::u,  {"u", "[m/s]"}},
    // </TODO>
    {solver_t::ix::rc, {"rc", "[kg kg-1]"}},
    {solver_t::ix::rr, {"rr", "[kg kg-1]"}},
    {solver_t::ix::nc, {"nc", "[kg-1]"}},
    {solver_t::ix::nr, {"nr", "[kg-1]"}}
  };
}
