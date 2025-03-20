#ifndef __IMPL_LOAD_BALANCER_HPP__
#define __IMPL_LOAD_BALANCER_HPP__

#include <load_balancing/iload_balancer.hpp>
#include <vector>

class UniformLoadBalancer final : public ILoadBalancer
{
public:
  explicit UniformLoadBalancer(uint32_t s);

protected:
  [[nodiscard]] double getRatio(uint64_t n, uint32_t rank) const noexcept final;
};

class HostImportantLoadBalancer final : public ILoadBalancer
{
public:
  HostImportantLoadBalancer(uint32_t s, double _alpha);

protected:
  [[nodiscard]] double getRatio(uint64_t n, uint32_t rank) const noexcept final;

private:
  double alpha;
};

class CustomLoadBalancer final : public ILoadBalancer
{
public:
  CustomLoadBalancer(uint32_t s, std::vector<double> _ratio);

protected:
  [[nodiscard]] double getRatio(uint64_t n, uint32_t rank) const noexcept final;

private:
  std::vector<double> ratios;
};

class BoundLoadBalancer final : public ILoadBalancer
{
public:
  BoundLoadBalancer(uint32_t s, uint64_t _n_max);

protected:
  [[nodiscard]] double getRatio(uint64_t n, uint32_t rank) const noexcept final;

private:
  uint64_t n_max;
};

#endif