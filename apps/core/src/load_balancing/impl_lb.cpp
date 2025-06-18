#include <load_balancing/impl_lb.hpp>
#include <stdexcept>
#include <utility>

UniformLoadBalancer::UniformLoadBalancer(uint32_t s) : ILoadBalancer(s) {};
[[nodiscard]] double UniformLoadBalancer::getRatio(uint64_t n, uint32_t rank) const noexcept
{
  (void)rank;
  (void)n;
  return 1. / static_cast<double>(size());
}

HostImportantLoadBalancer::HostImportantLoadBalancer(uint32_t s, double _alpha)
    : ILoadBalancer(s), alpha(_alpha)
{
  if (size() == 1)
  {
    alpha = 1.; // Overwrite user value to ensure correct calculation
  }
  else
  {
    if (alpha <= 0)
    {
      std::invalid_argument("Alpha should be strictly positive");
    }
  }
}

[[nodiscard]] double HostImportantLoadBalancer::getRatio(uint64_t n, uint32_t rank) const noexcept
{
  (void)n;
  const double ri = 1. / (alpha + size() - 1);

  return (rank == 0) ? alpha * ri : ri;
}

CustomLoadBalancer::CustomLoadBalancer(uint32_t s, std::vector<double> _ratio)
    : ILoadBalancer(s), ratios(std::move(_ratio))
{
  if (ratios.size() != static_cast<size_t>(size()))
  {
    throw std::invalid_argument("Ratio size should be equal to number of MPI nodes");
  }
}
[[nodiscard]] double CustomLoadBalancer::getRatio(uint64_t n, uint32_t rank) const noexcept
{
  (void)n;
  return ratios[rank];
}

BoundLoadBalancer::BoundLoadBalancer(uint32_t s, uint64_t _n_max) : ILoadBalancer(s), n_max(_n_max)
{
}

[[nodiscard]] double BoundLoadBalancer::getRatio(uint64_t n, uint32_t rank) const noexcept
{
  if (n_max > n)
  {
    return 0; // TODO Handle error
  }
  //just solve system where N nparticle and M n_max
  //N=M+(s-1)*alpha*N and M = beta*N


  double beta = static_cast<double>(n_max)/static_cast<double>(n);
  double alpha = (1. - beta) / static_cast<double>(size()-1);

  return (rank == 0) ? beta : alpha;
}
