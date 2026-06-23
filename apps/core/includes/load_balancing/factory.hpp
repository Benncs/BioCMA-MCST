#ifndef __CORE_LOAD_BALANCER_FACTORY_HPP__
#define __CORE_LOAD_BALANCER_FACTORY_HPP__

#include <load_balancing/iload_balancer.hpp>

std::unique_ptr<ILoadBalancer> lb_factory(uint32_t s);

#endif
