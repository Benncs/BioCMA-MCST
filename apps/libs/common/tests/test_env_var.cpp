#include <cassert>
#include <common/env_var.hpp>
#include <iostream>
#include <string>
#include <string_view>

template <typename T>
void
check(std::string_view name, T value, T vd)
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

void
check_read_env_optional(std::string_view name)
{

  assert(!Common::read_env<int>(name).has_value());

  assert(Common::set_local_env(name, 42));
  const auto value = Common::read_env<int>(name);
  assert(value.has_value());
  assert(*value == 42);

  assert(Common::set_local_env(name, std::string("not_a_number")));
  assert(!Common::read_env<int>(name).has_value());
  assert(Common::read_env_or<int>(name, 7) == 7);

#ifndef _WIN32
  unsetenv(name.data());
#endif
}

void
check_bool_spelling(std::string_view name)
{
  //works for true or 1
  for (const auto* spelling : { "true", "1" })
  {
    assert(Common::set_local_env(name, std::string(spelling)));
    assert(Common::read_env_or(name, false));
  }

  for (const auto* spelling : { "false", "0" })
  {
    assert(Common::set_local_env(name, std::string(spelling)));
    assert(!Common::read_env_or(name, true));
  }

  // Anything is false
  assert(Common::set_local_env(name, std::string("maybe")));
  assert(Common::read_env_or(name, true));
  assert(!Common::read_env_or(name, false));

#ifndef _WIN32
  unsetenv(name.data());
#endif
}

void
check_string_value(std::string_view name)
{
  assert(Common::set_local_env(name, std::string("ssavva")));
  assert(Common::read_env_or(name, std::string("fallback")) == "ssavva");

#ifndef _WIN32
  unsetenv(name.data());
  assert(Common::read_env_or(name, std::string("fallback")) == "fallback");
#endif
}

int
main()
{
  std::string test_var = "SSSAVAVA";
  check<int>(test_var, 1, 2);
  check<uint8_t>(test_var, 1, 2);
  check<uint64_t>(test_var, 1, 2);

  check<float>(test_var, 1., 2.);
  check<double>(test_var, 1., 2.);
  check<bool>(test_var, false, true);

  check_read_env_optional(test_var);
  check_bool_spelling(test_var);
  check_string_value(test_var);

  return 0;
}
