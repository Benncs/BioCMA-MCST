#include "load_balancing/impl_lb.hpp"
#include <cassert>
#include <iostream>
#include <load_balancing/iload_balancer.hpp>

void
test(ILoadBalancer* lb, uint64_t n, uint32_t n_rank)
{
  uint64_t cumsum = 0;
  for (uint32_t i = 0; i < n_rank; ++i)
  {

    auto ib = lb->balance(i, n);
    cumsum += ib;
    std::cerr << ib << std::endl;
  }
  std::cerr << cumsum << " " << n << std::endl;
  assert(lb->check(n));
  assert(cumsum == n);
}

#include <random>
int
main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(1, 10e6);

  for (auto n_rank = 1; n_rank < 10; ++n_rank)
  {
    for (auto _i = 0; _i < 10; ++_i)
    {
      const uint32_t n_particle = distrib(gen);
      UniformLoadBalancer uniform(n_rank);
      test(&uniform, n_particle, n_rank);

      HostImportantLoadBalancer host1(n_rank, 0.5);
      HostImportantLoadBalancer host2(n_rank, 3.1);

      test(&host1, n_particle, n_rank);
      test(&host2, n_particle, n_rank);
    }
  };
  auto n_rank = 10;
  CustomLoadBalancer custom(
      n_rank, { 0.2, 0.15, 0.1, 0.1, 0.2, 0.1, 0.05, 0.05, 0.025, 0.025 });
  BoundLoadBalancer bound(n_rank, 50);
  BoundLoadBalancer bound2(3, 4e6);
  const uint32_t n_particle = 10e6;

  UniformLoadBalancer uniform(6);
  test(&uniform, 10e3, 6);

  test(&custom, n_particle, n_rank);
  test(&bound, n_particle, n_rank);
  test(&bound2, 5e6, 3);
}
