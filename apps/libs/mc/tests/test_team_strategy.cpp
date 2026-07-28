#include "Kokkos_Macros.hpp"
#ifdef NDEBUG
#  undef NDEBUG
#endif

#include <Kokkos_Assert.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Printf.hpp>
#include <common/kokkos_getpolicy.hpp>

struct tag1
{
};
struct tag2
{
};

struct tag3
{
};

struct tag4_tile
{
};

struct tag_mlp
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

template <typename exec_space> struct f_assert
{

  using team_policy_t = Kokkos::TeamPolicy<exec_space>;
  using team_member_t = team_policy_t::member_type;
  using _float_t = float;
  using scratch_view_t = Kokkos::View<
      _float_t*,
      typename team_policy_t::execution_space::scratch_memory_space>;

  std::size_t n_tot;
  std::size_t m_n_per_team;

  explicit f_assert(std::size_t _n, std::size_t n_per_team)
      : n_tot(_n), m_n_per_team(n_per_team)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag_mlp _t,
             const team_member_t& team,
             std::size_t& reduce_val) const
  {
    (void)_t;
    const std::size_t N = m_n_per_team;
    const std::size_t p0 = team.league_rank() * N;

    const auto _ntot = n_tot;
    KOKKOS_ASSERT(p0 < _ntot)
    const auto upper_bound = ((p0 + N) >= n_tot) ? n_tot - p0 : N;

    KOKKOS_ASSERT(upper_bound > 0 && upper_bound <= n_tot);

    const std::size_t p = team.team_size();

    const std::size_t m = (upper_bound + p - 1) / p;

    std::size_t counter = 0;
    constexpr std::size_t U = 4;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 0, p),
        [&](const std::size_t tid, std::size_t& local)
        {
          for (std::size_t k = 0; k < m; k += U)
          {
            for (std::size_t u = 0; u < U; ++u)
            {
              const std::size_t idx = tid + (k + u) * p;

              if (idx >= upper_bound)
              {
                break;
              }

              const std::size_t flat_index = p0 + idx;

              KOKKOS_ASSERT(flat_index < _ntot);

              // Perform calculation here
              local += 1;
            }
          }
        },
        counter);

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() { reduce_val += counter; });
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag4_tile _t,
             const team_member_t& team,
             std::size_t& reduce_val) const
  {
    (void)_t;
    const std::size_t N = m_n_per_team;
    const std::size_t p0 = team.league_rank() * N;

    const auto _ntot = n_tot;
    KOKKOS_ASSERT(p0 < _ntot)
    const auto upper_bound = ((p0 + N) >= n_tot) ? n_tot - p0 : N;

    KOKKOS_ASSERT(upper_bound > 0 && upper_bound <= n_tot);

    const std::size_t p = team.team_size();

    const std::size_t m = (upper_bound + p - 1) / p;

    std::size_t counter = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 0, p),
        [&](const std::size_t tid, std::size_t& local)
        {
          for (std::size_t k = 0; k < m; ++k)
          {
            const std::size_t idx = tid + k * p; // stride p
            const std::size_t flat_index = p0 + idx;
            if (idx >= upper_bound)
            {
              break;
            }

            KOKKOS_ASSERT(flat_index < _ntot);

            // Perform calculation here
            local += 1;
          }
        },
        counter);

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() { reduce_val += counter; });
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag1 _t, const team_member_t& team, std::size_t& reduce_val) const
  {
    (void)_t;
    const auto npt = m_n_per_team;
    const std::size_t p0 = team.league_rank() * npt;
    const auto _ntot = n_tot;
    std::size_t local = 0;
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
        local);
    Kokkos::single(Kokkos::PerTeam(team), [&]() { reduce_val += local; });
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag2 _t, const team_member_t& team, std::size_t& reduce_val) const
  {
    (void)_t;
    const auto npt = m_n_per_team;
    const std::size_t p0 = team.league_rank() * npt;
    const auto _ntot = n_tot;
    std::size_t local = 0;
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
        local);

    Kokkos::single(Kokkos::PerTeam(team), [&]() { reduce_val += local; });
  }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(tag3 _t, const team_member_t& team, std::size_t& reduce_val) const
  {
    (void)_t;
    const auto npt = m_n_per_team;
    const std::size_t p0 = team.league_rank() * npt;
    const auto _ntot = n_tot;

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npt),
                         [=](int i)
                         {
                           const std::size_t p = p0 + i;
                           if (p >= _ntot)
                           {
                             return;
                           }
                         });
    std::size_t local = 0;
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
        local);
    Kokkos::single(Kokkos::PerTeam(team), [&]() { reduce_val += local; });
  }
};

template <typename tag>
void
test()
{
  const std::size_t n_tot = 1e9;

  const std::size_t n_per_team = cn_per_team();
  KOKKOS_ASSERT(n_per_team % 2 == 0);

  const std::size_t league_size = Common::c_league_size(n_tot, n_per_team);

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

template <typename tag>
void
test_smaller()
{
  const std::size_t n_tot = 52;
  const std::size_t n_per_team = 1024;

  const std::size_t league_size = Common::c_league_size(n_tot, n_per_team);
  // n<n_perteam should return 1
  KOKKOS_ASSERT(league_size == 1);

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
  test<tag4_tile>();
  test<tag_mlp>();
  test_smaller<tag1>();
  test_smaller<tag2>();
  test_smaller<tag3>();
  test_smaller<tag4_tile>();
  test_smaller<tag_mlp>();
  ;
  Kokkos::finalize();
}
