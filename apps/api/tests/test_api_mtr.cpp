#include <api/api.hpp>
#include <api/api_raw.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <simulation/mass_transfer.hpp>
#include <string>
#include <string_view>

#include "common_test.hpp"

#define INIT init_handle_raw(argc, argv);
#define INIT_CPP Api::SimulationInstance::init(argc, argv);
#define PARAM make_params(cx, ft, dt, np, nex, 0);

using MtrSetter = int (*)(Handle);

static constexpr std::array<MtrSetter, 3> mtr_setters
    = { set_mtr_auto, set_mtr_flowmap_turbulence, set_mtr_flowmap_kla };

void
mock_prepare_apply(std::string_view path, Handle handle)
{
  Param params = PARAM;
  auto rdir = tmp_dir; // tmp_dir comes from common_test
  CHECK(register_parameters(handle, &params));
  CHECK(register_result_path(handle, rdir.c_str()));
  CHECK(register_model_name(handle, "None"))
  CHECK(register_cma_path(handle, path.data()))
}

void
test_branch_null()
{
  // We use NULL to mimic C behavior
  for (const auto& set_mtr : mtr_setters)
  {
    CHECK_FALSE(set_mtr(NULL));

    const auto* msg = get_last_error();
    std::string_view v_msg = msg;
    assert(!v_msg.empty());
  }
}

void
test_set_mtr(int argc, char** argv)
{
  for (const auto& set_mtr : mtr_setters)
  {
    Handle handle = INIT;
    assert(handle != nullptr);
    CHECK(set_mtr(handle));
    delete_handle(&handle);
  }
}

void
test_set_mtr_override(int argc, char** argv)
{
  Handle handle = INIT;
  assert(handle != nullptr);
  for (const auto& set_mtr : mtr_setters)
  {
    CHECK(set_mtr(handle));
  }
  CHECK(set_mtr_auto(handle));
  delete_handle(&handle);
}

// Auto is the default, apply doesn´t need an explicit registration
void
test_apply_default(int argc, char** argv, std::string_view path)
{
  Handle handle = INIT;
  mock_prepare_apply(path, handle);
  CHECK(apply(handle, 0));
  delete_handle(&handle);
}

void
test_apply_set_mtr(int argc, char** argv, std::string_view path)
{
  for (const auto& set_mtr : mtr_setters)
  {
    Handle handle = INIT;
    mock_prepare_apply(path, handle);
    CHECK(set_mtr(handle));
    CHECK(apply(handle, 0));
    delete_handle(&handle);
  }
}

// FixedKla is not reachable from the C api because it carries a payload
void
test_set_mtr_cpp(int argc, char** argv)
{
  auto handle = *INIT_CPP;
  assert(handle != nullptr);
  assert(handle->set_mtr(Simulation::MassTransfer::Type::Auto{}));
  assert(handle->set_mtr(Simulation::MassTransfer::Type::FlowmapTurbulence{}));
  assert(handle->set_mtr(Simulation::MassTransfer::Type::FlowmapKla{}));
  assert(handle->set_mtr(Simulation::MassTransfer::Type::FixedKla{ { 0.2 } }));
}

int
main(int argc, char** argv)
{
  std::string cma_path = get_cma_path(argc, argv);

  test_branch_null();

  test_set_mtr(argc, argv);
  test_set_mtr_override(argc, argv);
  test_set_mtr_cpp(argc, argv);

  test_apply_default(argc, argv, cma_path);
  test_apply_set_mtr(argc, argv, cma_path);

  return 0;
}
