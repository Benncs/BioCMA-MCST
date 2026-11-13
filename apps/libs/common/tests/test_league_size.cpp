#include <cassert>
#include <common/kokkos_getpolicy.hpp>

// not power
static_assert(!Common::is_power_of_2(0));
static_assert(!Common::is_power_of_2(3));
static_assert(!Common::is_power_of_2(1025));

// power
static_assert(Common::is_power_of_2(1));
static_assert(Common::is_power_of_2(2));
static_assert(Common::is_power_of_2(1024));

void
test_fits_in_one_team()
{
  // n_tot <= n_per_team returns before the power of 2 assertion evrything will
  // return 1
  assert(Common::c_league_size(0, 1024) == 1);
  assert(Common::c_league_size(1, 1024) == 1);
  assert(Common::c_league_size(1024, 1024) == 1);
  assert(Common::c_league_size(10, 100) == 1);
}

void
test_ceil()
{
  // 1 team of 1024, 1 team of 1
  assert(Common::c_league_size(1025, 1024) == 2);
  assert(Common::c_league_size(2048, 1024) == 2);
  assert(Common::c_league_size(2049, 1024) == 3);
  assert(Common::c_league_size(100, 1) == 100);
  assert(Common::c_league_size(3, 2) == 2);
}

void
test_covers_every_particle()
{
  // The league must ensure that league_size * n_per_team>= n_tot
  constexpr std::size_t n_per_team = 256;
  for (std::size_t n_tot = 1; n_tot < 5000; ++n_tot)
  {
    const auto league = Common::c_league_size(n_tot, n_per_team);
    assert(league * n_per_team >= n_tot);
    assert((league - 1) * n_per_team < n_tot);
  }
}

int
main()
{
  test_fits_in_one_team();
  test_ceil();
  test_covers_every_particle();
}
