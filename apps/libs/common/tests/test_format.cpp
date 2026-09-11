#include <cassert>
#include <common/logger.hpp>

static constexpr std::string_view ref_1 = "Hello world!";

static constexpr std::string_view ref_2 = "Ssssava";

static constexpr std::string_view ref_3 = "2025 1 07 Sssaa";

[[maybe_unused]] static constexpr std::string_view ref_4 = "9.44 Sss";
#include <iostream>

int
main()
{
  // literal
  auto test_1 = IO::format("Hello", " world!");
  assert(test_1 == ref_1);
  std::cerr << "test_format_1" << std::endl;
  // hetero string type 1
  std::string p1 = "Ssss";
  std::string_view p2 = "ava";
  auto test_2 = IO::format(p1, p2);
  assert(test_2 == ref_2);
  std::cerr << "test_format_2" << std::endl;
  // hetero string type 2
  std::string_view p11 = "Ssss";
  std::string p22 = "ava";
  auto test_22 = IO::format(p11, p22);
  assert(test_22 == ref_2);
  std::cerr << "test_format_3" << std::endl;

  int y = 2025; // NOLINT
  uint64_t m = 1;
  constexpr std::string_view d = "07 ";

  auto test3 = IO::format(y, " ", m, " ", d, "Sssaa");
  std::cerr << "test_format_4: " << test3 << std::endl;
  assert(test3 == ref_3);

  // Floating point is difficult to test because of number of digits
  //
  //  float p = 9.44; // NOLINT
  //  std::string t = "Sss";
  //  auto test4 = IO::format(p, " ", t);
  //  std::cerr << "test_format_5 " << test4 << std::endl;
  //  assert(test4 == ref_4);

  // Compare with to_string instead of a literal
  float p = 9.44; // NOLINT
  std::string t = "Sss";
  auto test4 = IO::format(p, " ", t);
  std::cerr << "test_format_5: " << test4 << std::endl;
  assert(test4 == std::to_string(p) + " " + t);

  // No argument at all is an empty string, not a crash
  assert(IO::format().empty());

  // A single argument of each kind
  assert(IO::format("only") == "only");
  assert(IO::format(7) == "7");
  assert(IO::format(std::string("s")) == "s");

  // Signed values keep their sign
  assert(IO::format(-42) == "-42");
  assert(IO::format("n=", -1, "/", 2) == "n=-1/2");

  std::cerr << "test_format_6" << std::endl;
}
