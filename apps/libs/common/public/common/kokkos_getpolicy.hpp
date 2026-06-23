#ifndef __COMMON_KOKKOS_GET_POLICY_HPP__
#define __COMMON_KOKKOS_GET_POLICY_HPP__

#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <common/env_var.hpp>

namespace Common
{

  inline constexpr bool
  is_power_of_2(std::size_t x)
  {
    return x != 0 && (x & (x - 1)) == 0;
  }

  std::size_t c_league_size(std::size_t n_tot, std::size_t n_per_team) noexcept;

  template <typename Tag = void>
  Kokkos::TeamPolicy<ComputeSpace, Tag>
  get_policy_team(std::size_t league_size = 1)
  {
    const auto _league_size = read_env_or("BIOMC_LEAGUE_SIZE", league_size);
    return Kokkos::TeamPolicy<ComputeSpace, Tag>(_league_size, Kokkos::AUTO);
  }

  template <typename Tag = void>
  Kokkos::TeamPolicy<ComputeSpace, Tag>
  get_policy_team_from_npt(std::size_t n_tot, std::size_t n_per_team)
  {
    return Kokkos::TeamPolicy<ComputeSpace, Tag>(
        c_league_size(n_tot, n_per_team), Kokkos::AUTO(), Kokkos::AUTO());
  }
} // namespace Common
#endif
