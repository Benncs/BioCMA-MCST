#ifndef __ILOAD_BALANCER_HPP__
#define __ILOAD_BALANCER_HPP__

#include <cmath>
#include <cstdint>
#include <memory>

class ILoadBalancer
{
public:
  virtual ~ILoadBalancer() = default;
  explicit ILoadBalancer(uint32_t s);
  ILoadBalancer() = delete;
  ILoadBalancer(const ILoadBalancer&) = delete;
  ILoadBalancer(ILoadBalancer&&) = delete;
  ILoadBalancer& operator=(ILoadBalancer&&) = delete;
  ILoadBalancer& operator=(const ILoadBalancer&) = delete;

  [[nodiscard]] bool check(uint64_t n=0) const;

  uint64_t balance(uint32_t rank, uint64_t n);

protected:
  [[nodiscard]] virtual double getRatio(uint64_t n, uint32_t rank) const noexcept =0;

  [[nodiscard]] inline auto size() const noexcept
  {
    return _size;
  }

private:
  [[nodiscard]] uint64_t get_alloc(uint64_t n, uint32_t rank) const noexcept;
  uint32_t _size;
};

std::unique_ptr<ILoadBalancer> lb_factory();

#endif