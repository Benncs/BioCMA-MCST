#ifndef __TEST_COMMON_TEST_HPP__
#define __TEST_COMMON_TEST_HPP__
#include <cassert>
#include <filesystem>

#define CHECK(__stmt__) assert((__stmt__) == 0);
#define CHECK_FALSE(__stmt__) assert((__stmt__) != 0);
static const auto tmp_dir = std::filesystem::temp_directory_path();

// constexpr size_t n_rank = 1;
// constexpr size_t i_rank = 0;
constexpr size_t id = 123;
constexpr size_t nt = 4;

constexpr double cx = 1;
constexpr double ft = 1;
constexpr double dt = 0.;
constexpr size_t np = 40;
constexpr size_t nex = 0;

void test_init(int argc, char** argv);
void test_exec(int argc, char** argv, std::string_view path);
void test_exec_err(int argc, char** argv);

void test_apply(int argc, char** argv, std::string_view path);
void test_apply_err(int argc, char** argv);
void test_register_result_path(int argc, char** argv);
void test_register_parameters(int argc, char** argv);

inline std::string get_cma_path(int argc, char** argv)
{
  if (argc != 2)
  {
    assert(false && "Need cma path");
  }

  return argv[1];
}

#endif //__TEST_COMMON_TEST_HPP__