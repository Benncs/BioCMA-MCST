#include <common/env_var.hpp>
#include <load_balancing/iload_balancer.hpp>
#include <load_balancing/impl_lb.hpp>
#include <memory>

std::unique_ptr<ILoadBalancer>
lb_factory(uint32_t s)
{
  auto bounded = Common::read_env<uint32_t>("BIOMC_LBBOUND");
  if (bounded)
  {
    return std::make_unique<BoundLoadBalancer>(s, *bounded);
  }
  else
  {
    return std::make_unique<UniformLoadBalancer>(s);
  }
}
