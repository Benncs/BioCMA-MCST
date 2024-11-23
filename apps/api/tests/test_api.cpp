#include "common_test.hpp"
#include <api/api.hpp>
#include <optional>

#define INIT Handle::init(n_rank, i_rank, id, nt);

Core::UserControlParameters gparams(std::string_view path)
{
  bool serde = false;
  return {.biomass_initial_concentration = cx,
          .final_time = ft,
          .delta_time = dt,
          .number_particle = np,
          .n_thread = nt,
          .number_exported_result = nex,
          .recursive = false,
          .force_override = true,
          .serde = serde,
          .initialiser_path = "",
          .model_name = "None",
          .results_file_name = tmp_dir,
          .cma_case_path = std::string(path),
          .serde_file = std::nullopt};
}

void test_exec(std::string_view path)
{
}

void test_apply(std::string_view path)
{
}

void test_apply_from_param_serde(std::string_view path)
{
  // TODO
}

void test_exec_from_param(std::string_view path)
{
  auto handle = *INIT;
  assert(handle->register_parameters(gparams(path)));
  assert(handle->apply(false));
  assert(handle->exec());
}

void test_apply_from_param(std::string_view path)
{
  auto handle = *INIT;
  assert(handle->register_parameters(gparams(path)));
  assert(handle->apply(false));
}

void test_register_parameters(std::string_view path)
{
  auto handle = *INIT;
  assert(handle->register_parameters(gparams(path)));
}

int main(int argc, char** argv)
{
  std::string cma_path = get_cma_path(argc, argv);
  test_init();
  test_apply_err();
  test_exec_err();
  test_register_parameters(cma_path);

  test_apply_from_param(cma_path);
  test_apply(cma_path);
  test_apply_from_param_serde(cma_path);

  // FIXME test_exec_from_param(cma_path);
}

void test_exec_err()
{
  auto handle = *INIT;
  assert(!handle->exec());
}

void test_init()
{
  auto handle = INIT;

  assert(handle.has_value() && "init 1");

  auto& ptr = *handle;

  assert(ptr != nullptr && "Init 2");
}

void test_register_result_path()
{
  auto handle = *INIT;
  assert(handle->register_result_path(tmp_dir.c_str()) == true);
}

void test_apply_err()
{
  auto handle = *INIT
  {
    auto rc = handle->apply(false);
    assert(rc.invalid());
  }
  {
    auto rc = handle->apply(true);
    assert(rc.invalid());
  }
}
