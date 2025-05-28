#ifndef __COMMON_EXPORT_HPP__
#define __COMMON_EXPORT_HPP__

#include "Kokkos_Macros.hpp"
#include <biocma_cst_config.hpp>
#include <common/alg.hpp>
#include <common/execinfo.hpp>
#include <common/traits.hpp>
#include <string_view>

using ComputeSpace = Kokkos::DefaultExecutionSpace;
using HostSpace = Kokkos::HostSpace;


#ifdef ENABLE_KOKKOS_PROFILING
#  include <Kokkos_Profiling_ScopedRegion.hpp>
#  define PROFILE_SECTION(__label_section__)                                                       \
    Kokkos::Profiling::ScopedRegion region(__label_section__);
#else
#  define PROFILE_SECTION(__label_section__) ;
#endif

#define EIGEN_INDEX(__VALUE__) static_cast<int>(__VALUE__)

#ifndef NDEBUG
#  include <iostream>
#  include <source_location>
class Canary
{
public:
  Canary() = delete;
  Canary(const Canary&) = delete;
  Canary(Canary&&) = delete;
  Canary& operator=(const Canary&) = delete;
  Canary& operator=(Canary&&) = delete;

  explicit Canary(std::string_view lbl, std::source_location location) : _lbl(lbl)
  {
    // std::cout << "\033[1;31m" << location.function_name() << "\033[0m " << lbl << std::endl;
    std::cout << location.function_name() << ": " << lbl << std::endl;
  }

  ~Canary()
  {
    std::cout << "END " << _lbl << std::endl;
  }

private:
  std::string _lbl;
};
#  define MkCanary(x) Canary(x, std::source_location());
#else
#  define MkCanary(x)
#endif

inline Kokkos::TeamPolicy<ComputeSpace> get_policy_auto(std::size_t range)
{

  Kokkos::TeamPolicy<ComputeSpace> _policy;

  int recommended_team_size = _policy.team_size_recommended(
      KOKKOS_LAMBDA(const Kokkos::TeamPolicy<ComputeSpace>::member_type& team_handle) {
        std::size_t idx =
            team_handle.league_rank() * team_handle.team_size() + team_handle.team_rank();
        if (idx >= range)
        {
          return;
        }
      },
      Kokkos::ParallelForTag());
  int league_size = (static_cast<int>(range) + recommended_team_size - 1) / recommended_team_size;

  _policy = Kokkos::TeamPolicy<ComputeSpace>(league_size, recommended_team_size);

  return _policy;
}

#endif //__COMMON_EXPORT_HPP__