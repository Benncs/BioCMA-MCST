#ifndef __COMMON_KOKKOS_GET_POLICY_HPP__
#define __COMMON_KOKKOS_GET_POLICY_HPP__

#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <common/env_var.hpp>

namespace
{
  template <typename T> struct TagDetector
  {
    KOKKOS_FUNCTION void operator()(
        T _tag,
        const Kokkos::TeamPolicy<ComputeSpace>::member_type& team_handle) const
    {
      (void)team_handle;
      (void)_tag;
    }

    // void operator()(
    //     T _tag,
    //     const Kokkos::TeamPolicy<ComputeSpace>::member_type& team_handle)
    //     const
    // {
    //   (void)_tag;
    //   (void)team_handle;
    // }
    TagDetector() = default;
  };
} // namespace

namespace Common
{

  template <typename FunctorType, typename Tag = int>
  [[deprecated]] Kokkos::TeamPolicy<ComputeSpace, Tag>
  get_policy(const FunctorType& f, std::size_t range, bool reduce = false)
  {
    (void)reduce;
    Kokkos::TeamPolicy<Tag> _policy;

    const auto recommended_team_size =
        _policy.team_size_recommended(f, Kokkos::ParallelForTag{});

    const auto team_size =
        Common::read_env_or("BIOMC_TEAM_SIZE", recommended_team_size);

    const auto league_size =
        (static_cast<int>(range) + team_size - 1) / team_size;

    return Kokkos::TeamPolicy<ComputeSpace, Tag>(league_size, team_size);
  }

  template <typename Tag = int>
  [[deprecated]] Kokkos::TeamPolicy<ComputeSpace, Tag>
  get_policy(std::size_t range, bool reduce = false)
  {
    auto f = TagDetector<Tag>();
    return get_policy<decltype(f), Tag>(f, range, reduce);
  }

  template <typename Tag = void>
  Kokkos::TeamPolicy<ComputeSpace, Tag>
  get_policy_team(std::size_t league_size = 1)
  {
    const auto _league_size = read_env_or("BIOMC_LEAGUE_SIZE", league_size);
    return Kokkos::TeamPolicy<ComputeSpace, Tag>(_league_size, Kokkos::AUTO);
  }
} // namespace Common
#endif
