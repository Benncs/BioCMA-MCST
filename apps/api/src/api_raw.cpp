#include <api/api.hpp>
#include <api/api_raw.h>
#include <api/results.hpp>
#include <common/console.hpp>
#include <common/logger.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <span>
#include <sstream>
#include <string>
#include <utility>

#define CHECK_HANDLE_OR_RETURN                                                 \
  if (handle == nullptr)                                                       \
  {                                                                            \
    set_error_msg(handle_null_err_msg);                                        \
    return -1;                                                                 \
  }

namespace
{
  std::string last_error_msg;
  bool has_error = false;

  void
  set_error_msg(std::string_view e)
  {
    last_error_msg = e;
    has_error = true;
  }

  constexpr auto lambda_ok = [](auto) { return 0; };

  int
  set_and_log_error(Handle handle, auto e)
  {
    set_error_msg(e);
    if (const auto& log = handle->get_logger(); log != nullptr)
    {
      log->error(IO::format(" ", e));
    }
    return -3;
  }

  constexpr int ID_VERIF = 2025;
  [[maybe_unused]] constexpr int f_true = 1;
  [[maybe_unused]] constexpr int f_false = 0;

  constexpr std::string_view handle_null_err_msg = "API handle is invalid ";
} // namespace

/*FFI Feed descriptor*/

int
version_is_compatible(int major, int minor, int dev)
{
  return Api::version_is_compatible(major, minor, dev) ? 0 : -1;
}

// exposed to api
const char*
get_last_error()
{
  if (!has_error)
  {
    return "";
  }
  has_error = false;
  return last_error_msg.c_str();
}

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

FeedHandle
new_linear_feed_descriptor(double flow, double df, uint64_t input_position)
{
  Simulation::Feed::FeedDescriptor* fd        // NOLINT
      = new Simulation::Feed::FeedDescriptor; // NOLINT
  fd->input_position = input_position;
  fd->output_position = input_position;
  fd->flow = 0.;
  fd->extra = Simulation::Feed::Linear{ flow, df };

  return fd;
}

int
register_mixture_composition(Handle handle, char** names, int n_species)
{

  CHECK_HANDLE_OR_RETURN

  if (names == nullptr)
  {
    return -2;
  }

  std::vector<std::string> species_names;
  species_names.reserve(n_species);
  for (int i = 0; i < n_species; ++i)
  {
    std::string_view current_s_n = names[i]; // NOLINT
    species_names.emplace_back(current_s_n);
  }

  return handle->register_mixture_composition(species_names)
      .match([](auto) { return 0; },
             [&](auto e) { return set_and_log_error(handle, e); });
}

int
add_feed_descriptor(Handle handle, FeedHandle fd, int gas)
{
  CHECK_HANDLE_OR_RETURN

  if (fd == nullptr)
  {
    return -2;
  }

  const auto phase = gas != 0 ? Phase::Gas : Phase::Liquid;

  return handle->add_feed(*fd, phase)
      .match(lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
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
delete_constant_feed_descriptor(FeedHandle* fd)
{
  if (fd != nullptr && *fd != nullptr)
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
  bool f_reaction = (params.f_reaction != 0);

  auto p = Core::UserControlParameters::m_default();
  p.biomass_initial_concentration = params.biomass_initial_concentration;
  p.final_time = params.final_time, p.delta_time = params.delta_time;
  p.number_particle = params.number_particle, p.n_thread = params.n_thread;
  p.number_exported_result = params.number_exported_result;
  p.force_override = force_override;
  p.load_serde = load_serde;
  p.save_serde = save_serde;
  p.uniform_mc_init = uniform_mc_init;
  p.f_reaction = f_reaction;
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
  const auto f_uniform_init = f_false;
  const auto f_reaction = f_true;
  const auto n_thread = 1; // TODO Remove
  return { biomass_initial_concentration,
           final_time,
           delta_time,
           number_particle,
           n_thread,
           number_exported_result,
           f_false,
           f_false,
           save,
           f_uniform_init,
           f_reaction };
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
  CHECK_HANDLE_OR_RETURN

  handle->set_auto_mtr(); // FIXME

  return handle->apply(to_load != 0)
      .match(lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
}

Handle
init_handle_raw(int argc, char** argv)
{
  auto opt_handle = Api::SimulationInstance::init(argc, argv);
  if (opt_handle.has_value())
  {

    std::unique_ptr<Api::SimulationInstance> handle = std::move(*opt_handle);
    if (handle->get_exec_info().current_rank == 0)
    {
      auto logger = std::make_shared<IO::Console>();
      logger->toggle_all();
      handle->set_logger(logger);
    }

    return handle.release();
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
  CHECK_HANDLE_OR_RETURN;

  if (handle->get_id() != ID_VERIF)
  {
    set_error_msg("API handle version is not valid");
    return -1;
  }

  return handle->exec().match(
      lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
}

/*
    REGISTER
*/

int
register_result_path(Handle handle, const char* c)
{
  CHECK_HANDLE_OR_RETURN;

  if (c == nullptr)
  {
    return -2;
  }

  return (handle->register_result_path(c)) ? 0 : -1;
}

int
register_cma_path(Handle handle, const char* c)
{
  CHECK_HANDLE_OR_RETURN;

  if (c == nullptr)
  {
    return -2;
  }

  return handle->register_cma_path(c).match(
      lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
}

int
register_serde(Handle handle, const char* c)
{

  CHECK_HANDLE_OR_RETURN;

  if (c == nullptr)
  {
    return -2;
  }

  return (handle->register_serde(c)) ? 0 : -1;
}

int
register_model_name(Handle handle, const char* c)
{
  CHECK_HANDLE_OR_RETURN;

  if (c == nullptr)
  {
    return -2;
  }

  return handle->register_model_name(c).match(
      lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
}

int
register_initializer_path(Handle handle, const char* c)
{
  CHECK_HANDLE_OR_RETURN;

  if (c == nullptr)
  {
    return -2;
  }

  return handle->register_initialiser_file_path(c).match(
      lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
}

int
set_scalar_buffer(Handle handle,
                  uint64_t rows,
                  uint64_t cols,
                  double* liquid,
                  double* gas_ptr)
{

  CHECK_HANDLE_OR_RETURN;

  if (liquid == nullptr)
  {
    return -2;
  }

  if (rows == 0 || cols == 0)
  {
    return -1;
  }

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

    return handle
        ->register_scalar_initiazer(
            Core::ScalarFactory::FullCase(rows, std::move(liq), std::move(gas)))
        .match(lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
  }
  catch (const std::bad_alloc& e)
  {
    set_error_msg(e.what());
    return -1;
  }
  catch (const std::exception& e)
  {
    set_error_msg(e.what());
    return -1;
  }
}

int
register_parameters(Handle handle, Param* raw_params)
{
  CHECK_HANDLE_OR_RETURN;

  if (raw_params == nullptr)
  {
    return -2;
  }

  auto params = convert_c_wrap_to_param(*raw_params);
  return handle->register_parameters(std::move(params))
      .match(lambda_ok, [&](auto e) { return set_and_log_error(handle, e); });
}

// Parameters

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
  CHECK_HANDLE_OR_RETURN;
  return static_cast<int>(handle->get_exec_info().n_rank);
}
int
i_rank(Handle handle)
{
  CHECK_HANDLE_OR_RETURN;
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
