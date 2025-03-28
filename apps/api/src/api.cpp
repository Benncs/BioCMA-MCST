#include "api/results.hpp"
#include "common/common.hpp"
#include <Kokkos_Core.hpp>
#include <api/api.hpp>
#include <biocma_cst_config.hpp>
#include <common/execinfo.hpp>
#include <core/case_data.hpp>
#include <core/global_initaliser.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <memory>
#include <new>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
constexpr int ID_VERIF = 2025;

// #ifdef NO_MPI
// namespace WrapMPI
// {
//   bool is_initialized() noexcept
//   {
//     return false;
//   }
// } // namespace WrapMPI
// #else
// #  include "mpi_w/message_t.hpp"
// #  include <mpi_w/wrap_mpi.hpp>
// #endif

namespace
{
  bool check_required(const Core::UserControlParameters& params, bool to_load)
  {
    bool flag = true;

    flag = flag && params.final_time > 0;
    flag = flag && params.delta_time >= 0;
    flag = flag && !params.cma_case_path.empty();
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

} // namespace

namespace Api
{

  void finalise()
  {
    if (!Kokkos::is_finalized() && Kokkos::is_initialized())
    {
      Kokkos::finalize();
    }
  }

  [[nodiscard]] int SimulationInstance::get_id() const
  {
    return id;
  }

  template <typename T> std::vector<T> span_to_vec(std::span<T> rhs)
  {
    return std::vector<T>(rhs.begin(), rhs.end());
  }

  bool SimulationInstance::set_feed_constant_from_rvalue(double _f,
                                                         std::vector<double>&& _target,
                                                         std::vector<std::size_t>&& _position,
                                                         std::vector<std::size_t>&& _species,
                                                         bool gas,
                                                         bool fed_batch)
  {
    return set_feed_constant(_f, _target, _position, _species, gas, fed_batch);
  }

  bool SimulationInstance::set_feed_constant(double _flow,
                                             std::span<double> _concentration,
                                             std::span<std::size_t> _position,
                                             std::span<std::size_t> _species,
                                             bool gas,
                                             bool fed_batch)
  {
    auto target = span_to_vec(_concentration);
    auto position = span_to_vec(_position);
    auto species = span_to_vec(_species);
    auto fd = Simulation::Feed::FeedFactory::constant(
        _flow, std::move(target), std::move(position), std::move(species), !fed_batch);
    // negates fed_batch because constant accepts set_exit flag. fed_batch = !set_exit

    if (!feed.has_value())
    {
      feed = Simulation::Feed::SimulationFeed{std::nullopt, std::nullopt};
    }

    auto& current_descriptor = (gas) ? this->feed->gas : this->feed->liquid;

    // If not already created, new vector
    if (!current_descriptor.has_value())
    {
      current_descriptor = {fd};
    }
    else
    {
      // Else push
      current_descriptor->emplace_back(fd);
    };

    return true;
  }

  SimulationInstance::SimulationInstance(int argc, char** argv, std::optional<std::size_t> run_id)
      : id(ID_VERIF)
  {
    _data.exec_info = Core::runtime_init(argc, argv, run_id);
    std::atexit(Api::finalise);
  }

  SimulationInstance::~SimulationInstance()
  {
    _data = Core::CaseData(); // Explicity delete everything before
  }

 

  std::optional<std::unique_ptr<SimulationInstance>>
  SimulationInstance::init(int argc, char** argv, std::optional<std::size_t> run_id) noexcept
  {
    PROFILE_SECTION("Initialisation")
    auto ptr = std::unique_ptr<SimulationInstance>(new (std::nothrow)
                                                       SimulationInstance(argc, argv, run_id));
    if (ptr == nullptr)
    {
      return std::nullopt;
    }
    return ptr;
  }



  ApiResult SimulationInstance::exec() noexcept
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
      return ApiResult("Error apply first");
    }
  }

  ApiResult SimulationInstance::apply_load() noexcept
  {
    if (!check_required(this->params, true))
    {
      return ApiResult("Check params");
    }
    if (loaded)
    {
      return ApiResult("Not loaded");
    }
    if (auto opt_case = Core::load(this->_data.exec_info, std::move(this->params), this->feed))
    {
      this->_data = std::move(*opt_case);
      this->loaded = true;
      return ApiResult();
    }
    return ApiResult("Error loading");
  }

  ApiResult SimulationInstance::apply() noexcept
  {

    if (!check_required(this->params, false))
    {
      return ApiResult("Check params");
    }
    if (loaded)
    {
      return ApiResult("Not loaded");
    }
    if (!registered)
    {
      return ApiResult("Register first");
    }

    Core::GlobalInitialiser global_initializer(_data.exec_info, params);
    auto t = global_initializer.init_transitionner();
    global_initializer.init_feed(feed);

    auto __simulation = global_initializer.init_simulation(this->scalar_initializer_variant);
    if ((!t.has_value() && !__simulation.has_value()) || !global_initializer.check_init_terminate())
    {
      ApiResult("Error apply");
    }
    _data.params = global_initializer.get_parameters();
    _data.simulation = std::move(*__simulation);
    _data.transitioner = std::move(*t);
    applied = true;

    return ApiResult();
  }

  ApiResult SimulationInstance::register_scalar_initiazer(Core::ScalarFactory::ScalarVariant&& var)
  {
    this->scalar_initializer_variant = std::move(var);
    return ApiResult();
  }

  ApiResult SimulationInstance::apply(bool to_load) noexcept
  {
    if (to_load)
    {
      return apply_load();
    }

    return apply();
  }

  ApiResult SimulationInstance::register_parameters(Core::UserControlParameters&& _params) noexcept
  {

    params = std::move(_params);
    registered = true;
    return ApiResult();
  }

  bool SimulationInstance::register_result_path(std::string_view path)
  {

    // TODO Check path
    this->params.results_file_name = path;

    return true; // TODO
  }

  ApiResult SimulationInstance::register_initialiser_file_path(std::string_view path)
  {

    this->params.initialiser_path = path;
    return ApiResult(); // TODO
  }

  bool SimulationInstance::register_cma_path(std::string_view path, bool recursive)
  {
    // TODO Check path
    this->params.recursive = recursive;
    this->params.cma_case_path = path;
    return true; // TODO
  }

  ApiResult
  SimulationInstance::register_initial_condition(Core::ScalarFactory::ScalarVariant&& type)
  {
    scalar_initializer_variant = std::move(type);
    return ApiResult();
  }

  bool SimulationInstance::register_serde(std::string_view path)
  {
    this->params.serde_file = path;
    this->params.serde = true;
    return true;
  }

  ApiResult SimulationInstance::register_model_name(std::string_view path)
  {
    this->params.model_name = path;
    return ApiResult();
  }

  

} // namespace Api
