#include "common/console.hpp"
#include "common/logger.hpp"
#include "simulation/feed_descriptor.hpp"
#include <api/api.hpp>
#include <api/api_raw.h>
#include <api/results.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <utility>

constexpr int ID_VERIF = 2025;
[[maybe_unused]] constexpr int f_true = 1;
[[maybe_unused]] constexpr int f_false = 0;

/*FFI Feed descriptor*/

FeedHandle
new_constant_feed_descriptor(double flow, uint64_t input_position)
{
  Simulation::Feed::FeedDescriptor* fd        // NOLINT
      = new Simulation::Feed::FeedDescriptor; // NOLINT
  fd->input_position = input_position;
  fd->output_position = input_position;
  fd->flow = flow;
  fd->extra = Simulation::Feed::Constant{};

  return fd;
}

int
set_feed_descriptor(Handle handle, FeedHandle fd, int gas)
{
  const auto phase = gas != 0 ? Phase::Gas : Phase::Liquid;
  if (handle != nullptr)
  {
    if (fd != nullptr)
    {
      return handle->set_feed(*fd, phase) ? 0 : -3;
    }
    return -2;
  }
  return -1;
}

int
add_species(FeedHandle fh, double value, uint64_t i_species)
{
  int rc = -1;
  if (fh != nullptr)
  {
    fh->values.push_back({ value, i_species });
    rc = 0;
  }
  return rc;
}

int
set_fedbatch(FeedHandle fh)
{
  int rc = -1;
  if (fh != nullptr)
  {
    fh->output_position = std::nullopt;
    rc = 0;
  }
  return rc;
}

int
set_output_position(FeedHandle fh, uint64_t output_position)
{
  int rc = -1;
  if (fh != nullptr)
  {
    fh->output_position = output_position;
    rc = 0;
  }
  return rc;
}

int
delete_constant_feed_descriptor(FeedHandle** fd)
{
  if (fd != nullptr)
  {
    delete *fd; // NOLINT
    *fd = nullptr;
  }
  return 0;
}

/*FFI Parameters */

static Core::UserControlParameters
convert_c_wrap_to_param(const wrap_c_param_t& params)
{
  bool force_override = params.force_override != 0;
  bool load_serde = (params.load_serde != 0);
  bool save_serde = (params.save_serde != 0);
  bool uniform_mc_init = (params.uniform_particle_init != 0);

  auto p = Core::UserControlParameters::m_default();
  p.biomass_initial_concentration = params.biomass_initial_concentration;
  p.final_time = params.final_time, p.delta_time = params.delta_time;
  p.number_particle = params.number_particle, p.n_thread = params.n_thread;
  p.number_exported_result = params.number_exported_result;
  p.force_override = force_override;
  p.load_serde = load_serde;
  p.save_serde = save_serde;
  p.uniform_mc_init = uniform_mc_init;
  return p;
}

Param
make_params(double biomass_initial_concentration,
            double final_time,
            double delta_time,
            uint64_t number_particle,
            uint32_t number_exported_result,
            int save)
{

  return { biomass_initial_concentration,
           final_time,
           delta_time,
           number_particle,
           1,
           number_exported_result,
           f_false,
           f_false,
           save,
           f_false };
}

Param*
make_params_ptr(double biomass_initial_concentration,
                double final_time,
                double delta_time,
                uint64_t number_particle,
                uint32_t number_exported_result,
                int save)
{
  // NOLINTBEGIN(cppcoreguidelines-owning-memory)
  auto* const params = new Param(make_params(biomass_initial_concentration,
                                             final_time,
                                             delta_time,
                                             number_particle,
                                             number_exported_result,
                                             save));
  // NOLINTEND(cppcoreguidelines-owning-memory)
  return params;
}

void
delete_params(Param** params)
{
  if (params != nullptr)
  {
    delete *params; // NOLINT
    *params = nullptr;
  }
}

void
repr_user_param(const wrap_c_param_t* params, char** repr)
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

/*FFI API*/

// void finalize()
// {
//   Api::finalise();
// }

int
apply(Handle handle, int to_load)
{
  if (handle != nullptr)
  {
    handle->set_auto_mtr();                // FIXME
    auto rc = handle->apply(to_load != 0); // TODO HANDLE ERROR
    if (!rc)
    {
      handle->get_logger()->error(IO::format(" ", rc.get()));
    }
    return rc.to_c_ret_code();
  }
  return -1;
}

Handle
init_handle_raw(int argc, char** argv)
{
  auto opt_handle = Api::SimulationInstance::init(argc, argv);
  if (opt_handle.has_value())
  {
    auto logger = std::make_shared<IO::Console>();
    logger->toggle_all();
    (*opt_handle)->set_logger(logger);
    return opt_handle->release();
  }
  return nullptr;
}

void
delete_handle(Handle* handle)
{
  if (handle != nullptr)
  {

    delete *handle; // NOLINT
    *handle = nullptr;
  }
}

int
exec(Handle handle)
{
  if (handle != nullptr)
  {
    if (handle->get_id() == ID_VERIF)
    {
      auto rc = handle->exec();

      return rc ? 0 : -1;
    }
    return -2;
  }
  return -3;
}

/*
    REGISTER
*/

int
register_result_path(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_result_path(c)) ? 0 : -1;
  }
  return -1;
}

int
register_cma_path(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_cma_path(c)) ? 0 : -1;
  }
  return -1;
}

int
register_serde(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_serde(c)) ? 0 : -1;
  }
  return -1;
}

int
register_model_name(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_model_name(c)) ? 0 : -1;
  }
  return -1;
}

int
register_initializer_path(Handle handle, const char* c)
{
  if (handle != nullptr && c != nullptr)
  {
    return (handle->register_initialiser_file_path(c)) ? 0 : -1;
  }
  return -1;
}

int
set_scalar_buffer(Handle handle,
                  uint64_t rows,
                  uint64_t cols,
                  double* liquid,
                  double* gas_ptr)
{
  if (handle == nullptr || liquid == nullptr)
  {
    return -1;
  }

  if (rows == 0 || cols == 0)
  {
    return -1;
  }

  // uint64_t buffer_size=0;
  // if (__builtin_umull_overflow(rows, cols, &buffer_size)) {
  //     return -1;
  // }
  const auto buffer_size = rows * cols;

  try
  {
    std::span<double> liquid_span(liquid, buffer_size);
    std::vector<double> liq(liquid_span.begin(), liquid_span.end());

    std::optional<std::vector<double>> gas = std::nullopt;
    if (gas_ptr != nullptr)
    {
      std::span<double> gas_span(gas_ptr, buffer_size);
      gas = std::vector<double>(gas_span.begin(), gas_span.end());
    }

    bool success
        = handle
              ->register_scalar_initiazer(Core::ScalarFactory::FullCase(
                  rows, std::move(liq), std::move(gas)))
              .valid();

    return success ? 0 : -1;
  }
  catch (const std::bad_alloc& e)
  {
    return -1;
  }
  catch (const std::exception& e)
  {
    return -1;
  }
}

int
register_parameters(Handle handle, Param* raw_params)
{
  if (handle != nullptr && raw_params != nullptr)
  {
    auto params = convert_c_wrap_to_param(*raw_params);
    handle->register_parameters(std::move(params));
    return 0;
  }
  return -1;
}

int
set_feed_constant(Handle handle,
                  double flow,
                  double concentration,
                  size_t species,
                  size_t position,
                  int output_position,
                  int gas,
                  int fed_batch)
{
  if (handle != nullptr)
  {
    ApiResult res;

    const auto out_index = output_position == 0
                               ? std::nullopt
                               : std::make_optional(output_position);

    const auto phase = gas != 0 ? Phase::Gas : Phase::Liquid;
    // Negates fed_batch because expect set_output wich is !fed_batch

    const auto constant_feed = Simulation::Feed::FeedFactory::constant(
        flow, concentration, species, position, out_index, fed_batch != 0);
    res = handle->set_feed(constant_feed, phase);

    return res ? 0 : -1;
  }
  return -1;
}

// PArameters

void
show_user_param(const wrap_c_param_t* params)
{
  if (params != nullptr)
  {
    std::cout << convert_c_wrap_to_param(*params);
  }
}

int
n_rank(Handle handle)
{
  if (handle == nullptr)
  {
    return 0;
  }
  return static_cast<int>(handle->get_exec_info().n_rank);
}
int
i_rank(Handle handle)
{
  if (handle == nullptr)
  {
    return 0;
  }
  return static_cast<int>(handle->get_exec_info().current_rank);
}

void
get_model_list(char*** names, int* n_model)
{
  if (names == nullptr || n_model == nullptr)
  {
    return; // safety check
  }

  std::vector<std::string> models = Api::SimulationInstance::get_model_list();
  *n_model = static_cast<int>(models.size());

  // Allocate array of char* pointers
  *names = static_cast<char**>(std::malloc(sizeof(char*) * (*n_model)));
  if (*names == nullptr)
  {
    *n_model = 0;
    return; // allocation failed
  }

  for (int i = 0; i < *n_model; ++i)
  {
    const std::string& s = models[i];
    (*names)[i] = static_cast<char*>(std::malloc(s.size() + 1));
    if (!(*names)[i])
    {
      // Allocation failed: free previously allocated memory
      for (int j = 0; j < i; ++j)
      {
        std::free((*names)[j]);
      }
      std::free(*names);
      *names = nullptr;
      *n_model = 0;
      return;
    }
    std::strcpy((*names)[i], s.c_str()); // safe, we allocated enough space
  }
}

void
free_model_list(char** names, int n_model)
{
  if (names == nullptr)
  {
    return;
  }

  for (int i = 0; i < n_model; ++i)
  {
    std::free(names[i]);
  }
  std::free(names);
}