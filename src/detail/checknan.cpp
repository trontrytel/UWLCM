#pragma once

#ifdef NDEBUG
#define nancheck(arr, name) ((void)0)
#else
#define nancheck(arr, name) {nancheck_hlprs::nancheck_hlpr(arr, name);}
#endif

#ifdef NDEBUG
#define nancheck2(arrcheck, arrout, name) ((void)0)
#else
#define nancheck2(arrcheck, arrout, name) {nancheck_hlprs::nancheck2_hlpr(arrcheck, arrout, name);}
#endif

#ifdef NDEBUG
#define negcheck(arr, name) ((void)0)
#else
#define negcheck(arr, name) {nancheck_hlprs::negcheck_hlpr(arr, name);}
#endif

#ifdef NDEBUG
#define negcheck2(arrcheck, arrout, name) ((void)0)
#else
#define negcheck2(arrcheck, arrout, name) {nancheck_hlprs::negcheck2_hlpr(arrcheck, arrout, name);}
#endif

#ifdef NDEBUG
// actually not to zero, but to 1e-10 (we need rv>0 in libcloud and cond substepping numerical errors colud lead to rv<0 if we would set it here to 0)
#define negtozero(arr, name) {arr = where(arr <= 0., 1e-10, arr);}
#else
#define negtozero(arr, name) {nancheck_hlprs::negtozero_hlpr(arr, name);}
#endif

#ifdef NDEBUG
#define src_limiter(arr_val, arr_src, dt, name) {arr_src = where(arr_val + dt * arr_arc < 0., 0.99 * arr_val / dt, arr_src);}
#else
#define src_limiter(arr_val, arr_src, dt, name) {nancheck_hlprs::src_limiter_hlpr(arr_val, arr_src, dt, name);}
#endif

#ifndef NDEBUG
namespace nancheck_hlprs
{
  template<class arr_t>
  void nancheck_hlpr(const arr_t &arr, const std::string &name)
  {
    if(!std::isfinite(sum(arr)))
    {
      #pragma omp critical
      {
        std::cerr << "A not-finite number detected in: " << name << std::endl;
        std::cerr << arr;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void nancheck2_hlpr(const arr_t &arrcheck, const arr_t &arrout, const std::string &name)
  {
    if(!std::isfinite(sum(arrcheck)))
    {
      #pragma omp critical
      {
        std::cerr << "A not-finite number detected in: " << name << std::endl;
        std::cerr << arrcheck;
        std::cerr << arrout;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negcheck_hlpr(const arr_t &arr, const std::string &name)
  {
    if(min(arr) < 0.)
    {
      #pragma omp critical
      {
        std::cerr << "A negative number detected in: " << name << std::endl;
        std::cerr << arr;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negcheck2_hlpr(const arr_t &arrcheck, const arr_t &arrout,const std::string &name)
  {
    if(min(arrcheck) < 0.)
    {
      #pragma omp critical
      {
        std::cout << "A negative number detected in: " << name << std::endl;
        std::cout << arrcheck;
        std::cout << arrout;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negtozero_hlpr(arr_t arr, const std::string &name)
  {
    auto minval = min(arr);
    if(minval < 0.)
    {
      #pragma omp critical
      {
        std::cerr << "A negative number " << minval <<" detected in: " << name << std::endl;
        std::cerr << "CHEATING: turning negative values to small positive values" << std::endl;
      }
      arr = where(arr <= 0., 1e-10, arr);
    }
  }

  template<class arr_t, typename real_t>
  void src_limiter_hlpr(arr_t arr_val, arr_t arr_src, real_t dt, const std::string &name)
  {
    auto minval = min(arr_val + dt * arr_src);
    if(minval < 0.)
    {
      #pragma omp critical
      {
        std::cerr << "A negative number " << minval <<" is possible when applying sources to: " << name << std::endl;
        std::cerr << "CHEATING: turning the source terms resulting in negative values to " << name << "/dt" << std::endl;
      }
      arr_src = where(arr_val + dt * arr_src < 0., -0.99 * arr_val / dt , arr_src);
    }
  }
};
#endif
