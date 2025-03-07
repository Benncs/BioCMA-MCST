#include <load_balancing/impl_lb.hpp>
#include <stdexcept>
#include <utility>

UniformLoadBalancer::UniformLoadBalancer(uint32_t s) : ILoadBalancer(s) {};
[[nodiscard]] double UniformLoadBalancer::getRatio(uint32_t rank) const noexcept
{
  (void)rank;
  return 1. / static_cast<double>(size());
}

HostImportantLoadBalancer::HostImportantLoadBalancer(uint32_t s, double _alpha)
    : ILoadBalancer(s), alpha(_alpha)
{
  if (size() == 0)
  {
    alpha = 1.; // Overwrite user value to ensure correct calculation
  }
  else
  {
    if (size() - 1 <= alpha || alpha <= 0)
    {
      std::invalid_argument("Alpha has to respect inequality: size-1>alpha>0 ");
    }
  }
}

[[nodiscard]] double HostImportantLoadBalancer::getRatio(uint32_t rank) const noexcept
{
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
[[nodiscard]] double CustomLoadBalancer::getRatio(uint32_t rank) const noexcept
{
  return ratios[rank];
}
