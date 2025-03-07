#include <cassert>
#include <cstdint>
#include <load_balancing/iload_balancer.hpp>
ILoadBalancer::ILoadBalancer(uint32_t s) : _size(s) {};

uint64_t ILoadBalancer::get_alloc(uint64_t n, uint32_t rank) const noexcept
{
  return static_cast<uint64_t>(static_cast<double>(n) * getRatio(rank));
}

bool ILoadBalancer::check()const
{
  constexpr double tolerance = 1e-12;
  double cs = 0.;
  for (decltype(_size) i = 0; i < _size; ++i)
  {
    cs += getRatio(i);
  }
  return (std::abs(cs - 1.) < tolerance);
}

uint64_t ILoadBalancer::balance(uint32_t rank, uint64_t n)
{
  uint64_t allocated = 0;

  if (rank != 0)
  {
    allocated = get_alloc(n, rank);
  }
  else
  {
    uint64_t total_allocated = 0;
    for (decltype(_size) i = 0; i < _size; ++i)
    {
      allocated = get_alloc(n, i);
      total_allocated += allocated;
    }
    assert(n>=total_allocated);
    uint64_t diff = n - total_allocated;

    allocated = get_alloc(n, rank);
    allocated += diff;
  }

  return allocated;
}
