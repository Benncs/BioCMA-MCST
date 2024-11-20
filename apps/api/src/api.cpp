#include "biocma_cst_config.hpp"
#include "common/execinfo.hpp"
#include "core/case_data.hpp"
#include "core/global_initaliser.hpp"
#include "core/simulation_parameters.hpp"
#include "simulation/feed_descriptor.hpp"

#include <Kokkos_Core.hpp>
#include <api/api.hpp>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
constexpr int ID_VERIF = 2025;

#ifdef NO_MPI
namespace WrapMPI
{
  bool is_initialized() noexcept
  {
    return false;
  }
} // namespace WrapMPI
#else
#  include "mpi_w/message_t.hpp"
#endif

// namespace Api
// {

static bool check_required(const Core::UserControlParameters& params, bool to_load)
{
  bool flag = true;

  flag = flag && params.final_time > 0;
  flag = flag && params.delta_time >= 0;
  flag = flag && params.cma_case_path != "";
  flag = flag && ((params.serde && params.serde_file.has_value()) ||
                  (!params.serde && !params.serde_file.has_value()));

  if (!to_load)
  {
    flag = flag && params.biomass_initial_concentration !=
                       0; // Biomass could be 0 at start but OK to say that X>0
    flag = flag && params.number_particle > 0;
    // flag = flag && params.initialiser_path != "";
    // flag = flag && params.model_name != "";
  }

  return flag;
}

[[nodiscard]] int Handle::get_id() const
{
  return id;
}

template <typename T> std::vector<T> span_to_vec(std::span<T> rhs)
{
  return std::vector<T>(rhs.begin(), rhs.end());
}

bool Handle::set_feed_constant(double _f,
                               std::span<double> _target,
                               std::span<std::size_t> _position,
                               std::span<std::size_t> _species,
                               bool gas)
{
  auto target = span_to_vec(_target);
  auto position = span_to_vec(_position);
  auto species = span_to_vec(_species);
  auto fd = Simulation::Feed::FeedFactory::constant(
      _f, std::move(target), std::move(position), std::move(species));

  if (!feed.has_value())
  {
    feed = Simulation::Feed::SimulationFeed{std::nullopt, std::nullopt};
  }

  auto& fv = (gas) ? feed->gas : feed->liquid;

  if (!fv.has_value())
  {
    fv = {fd};
  }
  else
  {
    fv->emplace_back(fd);
  };

  return true;
}

Handle::Handle(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_process)
    : id(ID_VERIF)
{
  _data.exec_info = {.n_rank = n_rank,
                     .current_rank = current_rank,
                     .thread_per_process = thread_per_process,
                     .verbose = true,
                     .run_id = id};

  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    if (!WrapMPI::is_initialized())
    {
      throw std::invalid_argument("INIT MPI FIRST");
    }
  }

  if (!Kokkos::is_initialized())
  {
    Kokkos::initialize(
        Kokkos::InitializationSettings()
            .set_disable_warnings(false)
            .set_num_threads(static_cast<int32_t>(_data.exec_info.thread_per_process))
            .set_map_device_id_by("random"));
  }
}

std::optional<std::unique_ptr<Handle>> Handle::init(uint32_t n_rank,
                                                    uint32_t current_rank,
                                                    uint64_t id,
                                                    uint32_t thread_per_proces) noexcept
{

  auto ptr = std::unique_ptr<Handle>(new (std::nothrow)
                                         Handle(n_rank, current_rank, id, thread_per_proces));
  if (ptr == nullptr)
  {
    return std::nullopt;
  }
  return ptr;
}

Handle::~Handle()
{
  std::cout << "destructing from c++" << std::endl;
}

ApiResult Handle::exec() noexcept
{
  if (loaded || (registered && applied))
  {
    try
    {
      std::cout << "Running " << this->_data.exec_info.current_rank << "..." << std::endl;
      Core::exec(std::forward<Core::CaseData>(this->_data));
      return ApiResult();
    }
    catch (std::exception& e)
    {
      return ApiResult(e.what());
    }
  }
  else
  {
    return ApiResult();
  }
}

ApiResult Handle::apply(bool to_load) noexcept
{
  if (!check_required(this->params, to_load))
  {
    return ApiResult("Check params");
  }
  if (loaded)
  {
    return ApiResult("Not loaded");
  }

  if (to_load)
  {
    if (auto opt_case = Core::load(this->_data.exec_info, std::move(this->params), feed))
    {
      this->_data = std::move(*opt_case);
      this->loaded = true;
      return ApiResult();
    }
    return ApiResult("Error loading");
  }

  if (!registered)
  {
    return ApiResult("Register first");
  }

  Core::GlobalInitialiser gi(_data.exec_info, params);
  auto t = gi.init_transitionner();
  gi.init_feed(feed);

  auto __simulation = gi.init_simulation();
  if ((!t.has_value() && !__simulation.has_value()) || !gi.check_init_terminate())
  {
    ApiResult("Error apply");
  }
  _data.params = gi.get_parameters();
  _data.simulation = std::move(*__simulation);
  _data.transitioner = std::move(*t);
  applied = true;
  return ApiResult();
}

// void Handle::register_parameters()
// {
//   if (loaded)
//   {
//     return;
//   }
//   this->params = Core::UserControlParameters::m_default();
//   this->params.final_time = 10;
//   this->params.initialiser_path = "";
//   this->params.model_name = "None";
//   this->params.delta_time = 0.1;
//   this->params.number_particle = 100;
//   this->params.results_file_name = "test_api_new";
//   this->params.number_exported_result = 50;
//   this->params.cma_case_path = "/home/benjamin/Documents/code/cpp/BioCMA-MCST/cma_data/0d_mono/";
//   this->params.serde = false;
//   // mock_param.flow_files = {mock_param.user_params.cma_case_path};

//   registered = true;
// }

ApiResult Handle::register_parameters(Core::UserControlParameters&& _params) noexcept
{
  params = std::move(_params);
  registered = true;
  return ApiResult();
}

bool Handle::register_result_path(std::string_view path)
{
  // TODO Check path
  this->params.results_file_name = path;
  return true; // TODO
}

bool Handle::register_cma_path(std::string_view path, bool recursive)
{
  // TODO Check path
  this->params.recursive = recursive;
  this->params.cma_case_path = path;
  return true; // TODO
}

bool Handle::register_serde(std::string_view path)
{
  this->params.serde_file = path;
  this->params.serde = true;
  return true;
}

bool Handle::register_model_name(std::string_view path)
{
  this->params.model_name = path;
  return true;
}

// } //namespace Api
