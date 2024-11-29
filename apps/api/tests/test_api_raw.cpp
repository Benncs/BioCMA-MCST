#include <api/api_raw.h>
#include <cassert>
#include <cstring>
#include <filesystem>
#include <string_view>

#include "common_test.hpp"
#define INIT init_handle_raw(n_rank, i_rank, id, nt);
#define PARAM make_params(cx, ft, dt, np, nex);

void test_init()
{
  Handle handle = INIT assert(handle != nullptr);
  delete_handle(handle);
}

void test_delete_handle()
{
  Handle handle = INIT delete_handle(handle);
}

void mock_prepre_apply(std::string_view path, Handle handle)
{
  Param params = PARAM;
  auto rdir = std::string(tmp_dir)+"/";
  CHECK(register_parameters(handle, &params));
  CHECK(register_result_path(handle, rdir.c_str()))
  CHECK(register_model_name(handle, "None"))
  if (path.size() != 0)
  {
    CHECK(register_cma_path(handle, path.data()))
  }
}
void test_exec(std::string_view path)
{
  Handle handle = INIT mock_prepre_apply(path, handle);
  CHECK(apply(handle, 0));
  int result = exec(handle);
  assert(result == 0);
  delete_handle(handle);
}

void test_apply(std::string_view path)
{
  Handle handle = INIT mock_prepre_apply(path, handle);
  CHECK(apply(handle, 0));
  delete_handle(handle);
}

void test_apply_err()
{
  Handle handle = INIT assert(apply(handle, 0) != 0); // THIS SHOULD RETURN ERROR
  delete_handle(handle);
}

void test_exec_err()
{
  Handle handle = INIT assert(exec(handle) != 0); // THIS SHOULD RETURN ERROR
  delete_handle(handle);
}

// TODO TEST APPLY/EXEC WITH LOAD

void test_register_result_path()
{
  Handle handle = INIT int result = register_result_path(handle, "path");
  register_cma_path(handle, "path");
  assert(result == 0);
  delete_handle(handle);
}

void test_make_params()
{
  Param params = PARAM assert(params.biomass_initial_concentration == cx);
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
  Handle handle = INIT Param params = PARAM int result = register_parameters(handle, &params);
  assert(result == 0);
  delete_handle(handle);
}

void test_register_cma_path_recursive()
{
  Handle handle = INIT int result = register_cma_path_recursive(handle, "path");
  assert(result == 0);
  delete_handle(handle);
}

void test_register_cma_path()
{
  Handle handle = INIT int result = register_cma_path(handle, "path");
  assert(result == 0);
  delete_handle(handle);
}

void test_register_serde()
{
  Handle handle = INIT;
  int result = register_serde(handle, "serde");
  assert(result == 0);
  delete_handle(handle);
}

void test_register_model_name()
{
  Handle handle = INIT;
  int result = register_model_name(handle, "model");
  assert(result == 0);
  delete_handle(handle);
}

void test_branch_null()
{
  // We use NULL to mimic C behavior
  CHECK_FALSE(apply(NULL, 0));

  CHECK_FALSE(exec(NULL));
  Handle handle = INIT;

  // Manual fuzziing Alternatively check ptr validity, has to return error
  CHECK_FALSE(register_model_name(NULL, "valid"));
  CHECK_FALSE(register_model_name(handle, NULL));
  CHECK_FALSE(register_model_name(NULL, NULL));

  CHECK_FALSE(register_cma_path(NULL, "valid"));
  CHECK_FALSE(register_cma_path(handle, NULL));
  CHECK_FALSE(register_cma_path(NULL, NULL));

  CHECK_FALSE(register_result_path(NULL, "valid"));
  CHECK_FALSE(register_result_path(handle, NULL));
  CHECK_FALSE(register_result_path(NULL, NULL));

  CHECK_FALSE(register_cma_path_recursive(NULL, "valid"));
  CHECK_FALSE(register_cma_path_recursive(handle, NULL));
  CHECK_FALSE(register_cma_path_recursive(NULL, NULL));

  CHECK_FALSE(register_serde(NULL, "valid"));
  CHECK_FALSE(register_serde(handle, NULL));
  CHECK_FALSE(register_serde(NULL, NULL));

  CHECK_FALSE(register_parameters(NULL, NULL));
  Param params = PARAM;
  CHECK_FALSE(register_parameters(NULL, &params));

  delete_handle(handle);
}

int main(int argc, char** argv)
{

  std::string cma_path = get_cma_path(argc, argv);

  test_init();
  test_delete_handle();

  test_exec_err();
  test_apply_err();

  test_branch_null();

  test_make_params();
  test_register_parameters();
  test_register_result_path();
  test_register_cma_path_recursive();
  test_register_cma_path();
  test_register_serde();
  test_register_model_name();

  test_apply(cma_path);
  // FIXME test_exec(cma_path);
  test_exec(cma_path);
  return 0;
}
