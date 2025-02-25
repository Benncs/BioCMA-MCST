#include <biocma_cst_config.hpp>
#include <cassert>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/neighbors.hpp>
#include <cmt_common/cma_case.hpp>
#include <core/global_initaliser.hpp>
#include <core/scalar_factory.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <iostream>
#include <mc/mcinit.hpp>
#include <memory>
#include <optional>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <transitionner/transitioner_factory.hpp>
#include <utility>
#include <vector>
#include <wrap_init_model_selector.hpp>

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#endif

#ifndef NDEBUG
#  define VERBOSE_ERROR std::cerr << "Error LINE: " << __LINE__ << std::endl;
#else
#  define VERBOSE_ERROR
#endif

namespace
{
  const Core::ScalarFactory::ScalarVariant DefaultIntialiser = Core::ScalarFactory::Uniform({1.});

  const Core::ScalarFactory::ScalarVariant DefaultIntialiserTPF = Core::ScalarFactory::Local(
      {5., 0.}, {0}, std::vector<double>({0, 0.}), std::vector<size_t>({0}));

  void recur_path(std::string_view rootPath, Core::SimulationParameters& params);
} // namespace

namespace Core
{
  GlobalInitialiser::GlobalInitialiser(const ExecInfo& _info, UserControlParameters _user_params)
      : info(_info), user_params(std::move(_user_params)), is_host(info.current_rank == 0)
  {
    std::unique_ptr<CmaUtils::FlowMapTransitionner> transitioner = nullptr;

    params = SimulationParameters::init(user_params);
    f_init_gas_flow = info.current_rank == 0 && params.is_two_phase_flow; // NOLINT
  }

  template <typename T> using OptionalPtr = GlobalInitialiser::OptionalPtr<T>;

  CmtCommons::cma_exported_paths_t
  GlobalInitialiser::get_path_files(const std::string& cma_case_path)
  {
    CmtCommons::CMACaseInfo cma_case = CmtCommons::CMACaseInfoReader::load_case(
        cma_case_path + std::string(AutoGenerated::expected_cma_case_name));
    t_per_flowmap = cma_case.time_per_flowmap;
    return cma_case.paths;
  }

  void GlobalInitialiser::init_feed(std::optional<Simulation::Feed::SimulationFeed> feed)
  {
    validate_step(InitStep::Feed);
    this->feed = std::move(feed);
  }
  
  std::optional<bool>
  GlobalInitialiser::init_state(std::unique_ptr<CmaRead::FlowIterator>& flow_handle)
  {
    if (!check_steps(InitStep::FlowIterator))
    {
      return std::nullopt;
    }
    if (is_host)
    {
      if (!host_init_state(flow_handle).has_value())
      {
        return std::nullopt;
      }
    }

    mpi_broadcast();

    liquid_neighbors.set_row_major();
    validate_step(InitStep::InitState);

    return true;
  }

  OptionalPtr<CmaUtils::FlowMapTransitionner> GlobalInitialiser::init_transitionner()
  {

    auto handle = init_flow_iterator();
    if (handle.has_value())
    {
      return init_transitionner(std::move(*handle));
    }

    return std::nullopt;
  }

  OptionalPtr<Simulation::SimulationUnit>
  GlobalInitialiser::init_simulation(std::optional<Core::ScalarFactory::ScalarVariant> variant)
  {
    if (!check_steps(InitStep::Feed))
    {
      return std::nullopt;
    }

    auto scalar = variant.has_value() ? init_scalar(std::move(*variant)) : init_scalar();
    auto mc = init_monte_carlo();

    if (scalar.has_value() && mc.has_value())
    {
      return init_simulation(std::move(*mc), *scalar);
    }
    return std::nullopt;
  }

  OptionalPtr<MC::MonteCarloUnit> GlobalInitialiser::init_monte_carlo()
  {
    if (!check_steps(InitStep::FlowIterator, InitStep::InitState))
    {
      return std::nullopt;
    }

    const auto i_model = AutoGenerated::get_model_index_from_name(user_params.model_name);
    double total_mass = 0.;

    auto mc_unit = AutoGenerated::wrap_init_model_selector(
        i_model,
        info,
        user_params.number_particle,
        liquid_volume,
        CmaUtils::FlowMapTransitionner::get_neighbors_view(liquid_neighbors),
        user_params.biomass_initial_concentration,
        total_mass);

    double correct_total_mass = total_mass;
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
      correct_total_mass = 0;
#ifndef NO_MPI
      correct_total_mass = WrapMPI::all_reduce(total_mass);
#endif
    }

    MC::post_init_weight(mc_unit, user_params.biomass_initial_concentration, correct_total_mass);

    validate_step(InitStep::MC);

    return mc_unit;
  }

  OptionalPtr<CmaRead::FlowIterator> GlobalInitialiser::init_flow_iterator()
  {

    if (!is_host)
    {
      validate_step(InitStep::FlowIterator);
      return nullptr;
    }

    const auto flowmap_paths = get_path_files(user_params.cma_case_path);
    params.flow_files.clear();
    if (user_params.recursive)
    {

      recur_path(user_params.cma_case_path, params);
    }
    else
    {
      params.flow_files.emplace_back(user_params.cma_case_path);
    }

    if (params.flow_files.empty())
    {
      throw std::invalid_argument("Missing files path");
    }

    CmaRead::ReactorState const* state = nullptr;
    std::unique_ptr<CmaRead::FlowIterator> flow_handle;
    try
    {
      flow_handle = std::make_unique<CmaRead::FlowIterator>(params.flow_files, flowmap_paths);
      if (flow_handle == nullptr)
      {
        throw std::runtime_error("Flow map are not loaded");
      }
      std::cout << "Flowmap loaded: " << flow_handle->size() << std::endl;

      state = &flow_handle->get_unchecked(0);
      if (state == nullptr)
      {
        throw std::runtime_error("Reactor not correclty initialised");
      }
      std::cout << "Flowmap loaded with " << state->n_compartments << " compartments" << std::endl;
    }
    catch (const std::exception& e)
    {
      std::stringstream err;
      err << "Error while reading files\t:" << e.what();
      throw std::runtime_error(err.str());
    }
    if (flow_handle == nullptr)
    {
      return std::nullopt;
    }

    validate_step(InitStep::FlowIterator);
    return flow_handle;
  }

  OptionalPtr<CmaUtils::FlowMapTransitionner>
  GlobalInitialiser::init_transitionner(std::unique_ptr<CmaRead::FlowIterator>&& flow_handle)
  {
    if (!check_steps(InitStep::FlowIterator))
    {
      return std::nullopt;
    }
    init_state(flow_handle);
    // Calculate the total number of time steps
    f_init_gas_flow = info.current_rank == 0 && params.is_two_phase_flow;
    const auto n_t = static_cast<size_t>(user_params.final_time / params.d_t) + 1;
    // Transitioner handles flowmap transition between time step, flowmaps are
    // only located in host but transitioner handles cache and receiving for
    // workers
    auto transitioner = CmaUtils::get_transitioner(CmaUtils::FlowmapTransitionMethod::Discontinuous,
                                                   params.n_different_maps,
                                                   params.n_per_flowmap,
                                                   n_t,
                                                   std::move(flow_handle),
                                                   f_init_gas_flow);

    validate_step(InitStep::Transitioner);

    return transitioner;
  }

  OptionalPtr<Simulation::SimulationUnit>
  GlobalInitialiser::init_simulation(std::unique_ptr<MC::MonteCarloUnit> _unit,
                                     const Simulation::ScalarInitializer& scalar_init)
  {
    if (!check_steps(
            InitStep::FlowIterator, InitStep::Transitioner, InitStep::InitState, InitStep::Feed))
    {
      VERBOSE_ERROR
      return std::nullopt;
    }

    if (_unit == nullptr)
    {
      VERBOSE_ERROR
      return std::nullopt;
    }
    if (!Core::ScalarFactory::sanitize(scalar_init))
    {
      VERBOSE_ERROR
      return std::nullopt;
    }

    validate_step(InitStep::MC);
    validate_step(InitStep::Scalar);

    auto simulation = std::make_unique<Simulation::SimulationUnit>(
        std::move(_unit), scalar_init, std::move(feed));

    // FIXME
    if (AutoGenerated::get_model_index_from_name(user_params.model_name) == 4 &&
        simulation->getDimensions().n_species < 4)
    {
      throw std::runtime_error("Model must have 4 species concentrations");
    }
    validate_step(InitStep::SimulationUnit);

    return simulation;
  }

  std::optional<Simulation::ScalarInitializer>
  GlobalInitialiser::init_scalar(Core::ScalarFactory::ScalarVariant&& variant)
  {

    
    auto scalar_init =
        Core::ScalarFactory::scalar_factory(f_init_gas_flow, gas_volume, liquid_volume, variant);
    if (info.current_rank != 0) // FIXME
    {
      scalar_init.gas_flow = false;
    }
    validate_step(InitStep::Scalar);
    return scalar_init;
  }

  std::optional<Simulation::ScalarInitializer> GlobalInitialiser::init_scalar()
  {
    if (!check_steps(InitStep::FlowIterator, InitStep::InitState))
    {
      return std::nullopt;
    }

    // We init gas flow specific object only on the host and if user defined two
    // phase flow simulation
    f_init_gas_flow = info.current_rank == 0 && params.is_two_phase_flow;

    Core::ScalarFactory::ScalarVariant arg;
    if (user_params.initialiser_path.empty())
    {
      std::cerr << "WARNING: using Default Initialiser" << std::endl;
      arg = params.is_two_phase_flow ? DefaultIntialiserTPF : DefaultIntialiser;
    }
    else
    {
      arg = Core::ScalarFactory::File(params.n_compartments, user_params.initialiser_path);
    }

    Simulation::ScalarInitializer scalar_init =
        Core::ScalarFactory::scalar_factory(f_init_gas_flow, gas_volume, liquid_volume, arg);

    if (info.current_rank != 0) // FIXME
    {
      scalar_init.gas_flow = false;
    }

    validate_step(InitStep::Scalar);

    return scalar_init;
  }

  void GlobalInitialiser::set_initial_number_particle(uint64_t np) noexcept
  {
    params.number_particle = np;
  }

  bool GlobalInitialiser::check_init_terminate() const
  {
    int i = 0;
    for (auto&& step : validated_steps)
    {
      if (!step)
      {
        std::cerr << "Initialisation step " << i << " not completed!" << std::endl;
        return false;
      }
      ++i;
    }

    return true;
  }

  [[nodiscard]] SimulationParameters GlobalInitialiser::get_parameters() const
  {
    return params;
  }

  std::optional<bool>
  GlobalInitialiser::host_init_state(std::unique_ptr<CmaRead::FlowIterator>& flow_handle)
  {
    if (flow_handle == nullptr)
    {
      throw std::runtime_error("Reactor not correclty initialised");
    }
    const auto* fstate = &flow_handle->get_unchecked(0);
    if (fstate == nullptr)
    {
      throw std::runtime_error("Reactor not correclty initialised");
    }
    params.n_compartments = fstate->n_compartments;
    liquid_neighbors = fstate->liquid_flow.getViewNeighors();

    if (fstate->gas_flow.is_empty())
    {
      params.is_two_phase_flow = false;
    }

    worker_neighbor_data =
        std::vector<size_t>(liquid_neighbors.data().begin(), liquid_neighbors.data().end());

    params.n_different_maps = flow_handle->size();

    liquid_volume = fstate->liquidVolume;
    gas_volume = fstate->gasVolume;

    if (liquid_volume.empty())
    {
      return std::nullopt;
    }

    if (liquid_volume.size() > 1)
    {
      // When multiple compartments are present, it implies the existence of
      // internal hydrodynamic time scales. To account for this, the simulation's
      // explicit time step is calculated to approximate a CFL condition, with the
      // formula: time_step = min(residence_time) / 100. This approach ensures
      // that the fluid movement between two steps is accurately represented
      // without losing flow information.

      params.d_t = (user_params.delta_time == 0.) ? flow_handle->MinLiquidResidenceTime() / 100.
                                                  : user_params.delta_time;
    }
    else
    {
      // IF 0D reactor, no hydrodynamic time scale (internal)
      params.d_t = user_params.delta_time;
    }

    // Define the duration of each flowmap and compute steps per flowmap
    params.t_per_flow_map = t_per_flowmap;

    // TODO fix float/integer division
    const auto n_per_flowmap =
        (t_per_flowmap == 0 || params.n_different_maps == 1)
            ? 1
            : static_cast<size_t>(t_per_flowmap / static_cast<double>(params.d_t)) + 1;

    params.n_per_flowmap = n_per_flowmap;

    return true;
  }

  void GlobalInitialiser::mpi_broadcast() // NOLINT
  {
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
#ifndef NO_MPI
      // Check once if we can broadcast integer, if not quit
      if (WrapMPI::broadcast(params.n_different_maps, 0) != 0)
      {
        WrapMPI::critical_error();
      }
      // broadcasting vector (malloc) is handled by WrapMPI so don't need to resize
      // before
      WrapMPI::broadcast(params.d_t, 0);
      WrapMPI::broadcast(params.n_per_flowmap, 0);
      WrapMPI::broadcast(params.n_compartments, 0);
      WrapMPI::broadcast(params.is_two_phase_flow, 0);
      WrapMPI::broadcast(liquid_volume, 0, info.current_rank);
      WrapMPI::broadcast(gas_volume, 0, info.current_rank);
      WrapMPI::broadcast(worker_neighbor_data, 0, info.current_rank);
      // WrapMPI::bcst_iterator(_flow_handle, info.current_rank);
      if (!is_host)
      {
        // Neighbors data is a matrix with nrow the number of compartment (known by
        // host and workers) n_col is not directly known but the matrix size is so
        // we calculate column this way
        assert(worker_neighbor_data.size() > params.n_compartments);
        const size_t n_col = worker_neighbor_data.size() / params.n_compartments;
        assert(worker_neighbor_data.size() % params.n_compartments == 0);

        liquid_neighbors = CmaRead::Neighbors::Neighbors_const_view_t(
            worker_neighbor_data, params.n_compartments, n_col);
      }
      if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
      {
        WrapMPI::barrier(); // This barrier is probably useless
      }

#endif
    }
  }

} // namespace Core

namespace
{

  void recur_path(std::string_view rootPath, Core::SimulationParameters& params)
  {
    size_t count = 1;
    std::string dirName = "i_" + std::to_string(count) + "/";
    std::filesystem::path dirPath = std::string(rootPath) + dirName;
    while (std::filesystem::exists(dirPath) && std::filesystem::is_directory(dirPath))
    {
      ++count;
      params.flow_files.push_back(dirPath.string());
      dirName = "i_" + std::to_string(count) + "/";
      dirPath = std::string(rootPath) + dirName;
    }
  }
} // namespace
