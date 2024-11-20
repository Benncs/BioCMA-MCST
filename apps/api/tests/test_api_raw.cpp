#include <api/api_raw.h>
#include <cassert>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string_view>

#define CHECK(__stmt__) assert((__stmt__) == 0);

static const auto tmp_dir = std::filesystem::temp_directory_path();

constexpr size_t n_rank = 1;
constexpr size_t i_rank = 0;
constexpr size_t id = 123;
constexpr size_t nt = 4;

constexpr double cx =1 ;
constexpr double ft =1 ;
constexpr double dt =0.;
constexpr size_t np =40;
constexpr size_t nex =0;
#define INIT init_handle_raw(n_rank, i_rank, id, nt);
#define PARAM make_params(cx,ft, dt, np, nex);

void test_init_handle_raw()
{
  Handle* handle = INIT
  assert(handle != nullptr);
  delete_handle(handle);
}

void test_delete_handle()
{
  Handle* handle = INIT
  delete_handle(handle);
}

void mock_prepre_apply(std::string_view path, Handle* handle)
{
  Param params = PARAM
  CHECK(register_parameters(handle, &params));
  CHECK(register_result_path(handle, tmp_dir.c_str()))
  CHECK(register_model_name(handle, "None"))
  CHECK(register_cma_path(handle, path.data()))
}
void test_exec(std::string_view path)
{
  Handle* handle = INIT
  mock_prepre_apply(path, handle);
  CHECK(apply(handle, 0));
  int result = exec(handle);
  assert(result == 0);
  delete_handle(handle);
}

void test_apply(std::string_view path)
{
  Handle* handle = INIT
  mock_prepre_apply(path, handle);
  CHECK(apply(handle, 0));
  delete_handle(handle);
}

//TODO TEST APPLY/EXEC WITH LOAD

void test_register_result_path()
{
  Handle* handle = INIT
  int result = register_result_path(handle, "path");
  register_cma_path(handle, "path");
  assert(result == 0); // Assuming 0 is a valid return value
  delete_handle(handle);
}

void test_make_params()
{
  Param params = PARAM
  assert(params.biomass_initial_concentration ==cx);
  assert(params.final_time == ft);
  assert(params.delta_time == dt);
  assert(params.number_particle == np);
  assert(params.n_thread == 1); // Default value
  assert(params.number_exported_result == nex);
  assert(params.recursive == 0);      // Default value
  assert(params.force_override == 0); // Default value
  assert(params.serde == 0);          // Default value
}

void test_register_parameters()
{
  Handle* handle = INIT
  Param params = PARAM
  int result = register_parameters(handle, &params);
  assert(result == 0); // Assuming 0 is a valid return value
  delete_handle(handle);
}

void test_register_cma_path_recursive()
{
  Handle* handle = INIT
  int result = register_cma_path_recursive(handle, "path");
  assert(result == 0); // Assuming 0 is a valid return value
  delete_handle(handle);
}

void test_register_cma_path()
{
  Handle* handle = INIT
  int result = register_cma_path(handle, "path");
  assert(result == 0); // Assuming 0 is a valid return value
  delete_handle(handle);
}

void test_register_serde()
{
  Handle* handle = INIT
  int result = register_serde(handle, "serde");
  assert(result == 0); // Assuming 0 is a valid return value
  delete_handle(handle);
}

void test_register_model_name()
{
  Handle* handle = INIT
  int result = register_model_name(handle, "model");
  assert(result == 0); // Assuming 0 is a valid return value
  delete_handle(handle);
}

int main(int argc, char** argv)
{

  if (argc != 2)
  {
    assert(false && "Need cma path");
  }

  std::string cma_path = argv[1]; //NOLINT

  test_init_handle_raw();
  test_delete_handle();

  test_make_params();
  test_register_parameters();
  test_register_result_path();
  test_register_cma_path_recursive();
  test_register_cma_path();
  test_register_serde();
  test_register_model_name();

  test_apply(cma_path);
  test_exec(cma_path);

  std::cout << "All tests passed!" << std::endl;
  return 0;
}
