#include <Kokkos_Assert.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Printf.hpp>

struct tag1
{
};
struct tag2
{
};

struct tag3
{
};

constexpr std::size_t
cn_per_team()
{
  const std::size_t l1_bytes = 32 * 1024;
  const int n_arrays = 3;
  const std::size_t elem_bytes = sizeof(double);

  const int m_n_per_team = l1_bytes / (n_arrays * elem_bytes);
  return std::pow(2, std::ceil(std::log(m_n_per_team) / std::log(2)));
}

std::size_t
cleague_size(std::size_t n_tot, std::size_t n_per_team)
{
  return (n_tot + n_per_team - 1) / n_per_team;
}

template <typename exec_space> struct f_assert
{

  using team_policy_t = Kokkos::TeamPolicy<exec_space>;
  using team_member_t = team_policy_t::member_type;

  using scratch_view_t = Kokkos::View<
      float_t*,
      typename team_policy_t::execution_space::scratch_memory_space>;
  using float_t = float;

  std::size_t n_tot;
  std::size_t m_n_per_team;

  explicit f_assert(std::size_t _n, std::size_t n_per_team)
      : n_tot(_n), m_n_per_team(n_per_team)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag1 _t, const team_member_t& team, std::size_t& reduce_val) const
  {
    (void)_t;
    const auto npt = m_n_per_team;
    const std::size_t p0 = team.league_rank() * npt;
    const auto _ntot = n_tot;
    std::size_t local_val = 0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, npt),
        [=](int i, std::size_t& lv)
        {
          const std::size_t p = p0 + i;
          if (p >= _ntot)
          {
            return;
          }
          lv += 1;
        },
        local_val);
    reduce_val += local_val;
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag2 _t, const team_member_t& team, std::size_t& reduce_val) const
  {
    (void)_t;
    const auto npt = m_n_per_team;
    const std::size_t p0 = team.league_rank() * npt;
    const auto _ntot = n_tot;
    std::size_t local_val = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, npt),
        [=](int i, std::size_t& lv)
        {
          const std::size_t p = p0 + i;
          if (p >= _ntot)
          {
            return;
          }
          lv += 1;
        },
        local_val);
    reduce_val += local_val;
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag3 _t, const team_member_t& team, std::size_t& reduce_val) const
  {
    (void)_t;
    const auto npt = m_n_per_team;
    const std::size_t p0 = team.league_rank() * npt;
    const auto _ntot = n_tot;
    std::size_t local_val = 0;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npt),
                         [=](int i)
                         {
                           const std::size_t p = p0 + i;
                           if (p >= _ntot)
                           {
                             return;
                           }
                         });
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, npt),
        [=](int i, std::size_t& lv)
        {
          const std::size_t p = p0 + i;
          if (p >= _ntot)
          {
            return;
          }
          lv += 1;
        },
        local_val);
    reduce_val += local_val;
  }
};

template <typename tag>
void
test()
{
  const std::size_t n_tot = 1e9;

  const std::size_t n_per_team = cn_per_team();
  KOKKOS_ASSERT(n_per_team % 2 == 0);

  const std::size_t league_size = cleague_size(n_tot, n_per_team);

  auto policy
      = Kokkos::TeamPolicy<tag>(league_size, Kokkos::AUTO(), Kokkos::AUTO());

  Kokkos::printf("N per team: %ld\r\n\tTeam "
                 "policy:\r\n\tleague_size=%d\r\n\tteam_size=%d\r\n",
                 n_per_team,
                 league_size,
                 policy.team_size());

  KOKKOS_ASSERT(league_size * n_per_team * policy.team_size() >= n_tot);
  std::size_t count = 0;
  Kokkos::parallel_reduce(
      policy,
      f_assert<Kokkos::DefaultExecutionSpace>(n_tot, n_per_team),
      count);

  Kokkos::printf("n=%ld\tcount=%ld\r\n", n_tot, count);

  KOKKOS_ASSERT(count == n_tot);
}

int
main()
{
  Kokkos::initialize();

  test<tag1>();
  test<tag2>();
  test<tag3>();

  Kokkos::finalize();
}