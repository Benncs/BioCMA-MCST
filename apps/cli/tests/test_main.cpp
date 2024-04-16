#include <cassert>
#include <common/common.hpp>

int main()
{
  test_common();

  static_assert(2 == (1 + 1));
  return 0;
}