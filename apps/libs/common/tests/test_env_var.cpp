#include <cassert>
#include <common/env_var.hpp>
#include <iostream>
#include <string_view>

template <typename T> void check(std::string_view name, T value, T vd)
{

  T fallback_result = Common::read_env_or(name, vd);
  assert(fallback_result == vd);
  bool set_result = Common::set_local_env(name, value);
  assert(set_result && "Failed to set environment variable");
  T result = Common::read_env_or(name, vd);
  std::cout << result << " " << vd << std::endl;
  assert(result == value);
#ifndef _WIN32
  unsetenv(name.data()); // This is to ensure there 's no side effect
#endif
}

int main()
{
  std::string test_var = "MY_TEST_VAR";
  check<int>(test_var, 1, 2);
  check<uint8_t>(test_var, 1, 2);
  check<uint64_t>(test_var, 1, 2);

  check<float>(test_var, 1., 2.);
  check<double>(test_var, 1., 2.);
  check<bool>(test_var, false, true);

  return 0;
}
