#include "core/scalar_factory.hpp"
#include <api/api.hpp>
#include <api/api_raw.h>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <optional>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <utility>
constexpr int ID_VERIF = 2025;

int apply(Handle handle, int to_load)
{
  if (handle != nullptr)
  {
    auto rc = handle->apply(to_load != 0); // TODO HANDLE ERROR
    return rc.to_c_ret_code();
  }
  return -1;
}

Handle init_handle_raw(int n_rank, int current_rank, uint64_t id, uint32_t thread_per_process)
{
  auto opt_handle = Api::SimulationInstance::init(n_rank, current_rank, id, thread_per_process);
  if (opt_handle.has_value())
  {
    return opt_handle->release();
  }
  return nullptr;
}

Handle init_handle_shared(uint64_t id, uint32_t thread_per_process)
{
  return init_handle_raw(1, 0, id, thread_per_process);
}


void delete_handle(Handle* handle)
{

  if (handle != nullptr)
  {
    delete *handle; // NOLINT
    *handle = nullptr;
  }
}

int exec(Handle handle)
{
  if (handle != nullptr)
  {
    if (handle->get_id() == ID_VERIF)
    {
      auto rc = handle->exec();
      std::cerr << "EXEC: " << rc.get() << std::endl;

      return rc ? 0 : -1;
    }
    return -2;
  }
  return -3;
}

/*
    REGISTER
*/

int register_result_path(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_result_path(c)) ? 0 : -1;
  }
  return -1;
}

int register_cma_path(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_cma_path(c, false)) ? 0 : -1;
  }
  return -1;
}

int register_cma_path_recursive(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {

    return (handle->register_cma_path(c, true)) ? 0 : -1;
  }
  return -1;
}

int register_serde(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_serde(c)) ? 0 : -1;
  }
  return -1;
}

int register_model_name(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_model_name(c)) ? 0 : -1;
  }
  return -1;
}

int register_initializer_path(Handle handle, const char* c)
{

  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_initialiser_file_path(c)) ? 0 : -1;
  }
  return -1;
}

Core::UserControlParameters convert_c_wrap_to_param(const wrap_c_param_t& params)
{
  bool recursive = params.recursive != 0;
  bool serde = params.serde != 0;
  bool force_override = params.force_override != 0;

  return {.biomass_initial_concentration = params.biomass_initial_concentration,
          .final_time = params.final_time,
          .delta_time = params.delta_time,
          .number_particle = params.number_particle,
          .n_thread = params.n_thread,
          .number_exported_result = params.number_exported_result,
          .recursive = force_override,
          .force_override = recursive,
          .serde = serde,
          .initialiser_path = "",
          .model_name = "",
          .results_file_name = "",
          .cma_case_path = "",
          .serde_file = std::nullopt};
}

Param make_params(double biomass_initial_concentration,
                  double final_time,
                  double delta_time,
                  uint64_t number_particle,
                  uint32_t number_exported_result)
{
  return {biomass_initial_concentration,
          final_time,
          delta_time,
          number_particle,
          1,
          number_exported_result,
          0,
          0,
          0};
}

int register_parameters(Handle handle, Param* raw_params)
{

  if (handle != nullptr && raw_params != nullptr)
  {
    auto params = convert_c_wrap_to_param(*raw_params);
    handle->register_parameters(std::move(params));
    return 0;
  }
  return -1;
}

int set_feed_constant(Handle handle,
                      double _f,
                      size_t n_species,
                      double* _target,
                      size_t* _species,
                      size_t n_position,
                      size_t* _position,
                      int gas,
                      int fed_batch)
{
  if ((handle != nullptr) && (_target != nullptr) && (_species != nullptr) &&
      (_position != nullptr))
  {
    auto span_target = std::span<double>(_target, n_species);
    auto span_species = std::span<std::size_t>(_species, n_species);
    auto span_pos = std::span<std::size_t>(_position, n_position);
    handle->set_feed_constant(_f, span_target, span_pos, span_species, gas != 0, fed_batch != 0);
    return 0;
  }

  return -1;
}
void show_user_param(const wrap_c_param_t* params)
{
  if (params != nullptr)
  {
    std::cout << convert_c_wrap_to_param(*params);
  }
}

void repr_user_param(const wrap_c_param_t* params, char** repr)
{
  std::stringstream ss;
  ss << convert_c_wrap_to_param(*params);

  // Allocate memory for the string
  std::string str = ss.str();
  // clang-format off
  *repr = static_cast<char*>(malloc((str.size() + 1) * sizeof(char))); // NOLINT +1 for null terminator
  // clang-format on
  if (*repr != nullptr)
  {
    strcpy(*repr, str.c_str()); // Copy the string to the allocated memory
  }
};
