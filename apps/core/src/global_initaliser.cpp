#include <algorithm>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <common/logger.hpp>
#include <core/global_initaliser.hpp>
#include <core/scalar_factory.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdio>
#include <dataexporter/data_exporter.hpp>
#include <load_balancing/iload_balancer.hpp>
#include <load_balancing/impl_lb.hpp>
#include <mc/mcinit.hpp>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/mass_transfer.hpp>
#include <simulation/simulation.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
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
  const Core::ScalarFactory::ScalarVariant DefaultIntialiser =
      Core::ScalarFactory::Uniform({1.});

  const Core::ScalarFactory::ScalarVariant DefaultIntialiserTPF =
      Core::ScalarFactory::Local({0., 0.},
                                 {0},
                                 std::vector<double>({0, 0.}),
                                 std::vector<size_t>({0}));

  // TODO Move elsewhere
  double get_time_step(double user_deta_time,
                       const CmaUtils::TransitionnerPtrType& iterator)
  {

    // internal hydrodynamic time scales. To account for this, the simulation's
    // explicit time step is calculated to approximate a CFL condition, with the
    // formula: time_step = min(residence_time) / 100. This approach ensures
    // that the fluid movement between two steps is accurately represented
    // without losing flow information.
    double delta_time = user_deta_time;

    if (delta_time <= 0)
    {

      const std::size_t n_states = iterator->size();
      double min_residence_time = std::numeric_limits<double>::max();
      for (std::size_t i_state = 0; i_state < n_states; ++i_state)
      {
        // Naive min-research, as COO is not sorted we can´t do better without
        // convert COO to CSR/CSC or sort
        // + function exececuted only once
        const auto state = iterator->get_at(i_state);
        const auto liquid = state->get_liquid();
        const auto coo_matrix = liquid->transition();
        const auto liquid_volumes = liquid->volume();
        const auto data = coo_matrix->values();
        const auto rows = coo_matrix->row_indices();
        const auto cols = coo_matrix->col_indices();
        for (std::size_t k = 0; k < data.size(); ++k)
        {
          // Only way to find diagonal with not sorted COO
          if (rows[k] == cols[k])
          {
            const auto volume = liquid_volumes[rows[k]];
            if (volume != 0.)
            {
              double residence_time =
                  -data[k] / volume; // minus cause transition has negative
                                     // diagonal (leaving flow)
              min_residence_time = std::min(residence_time, min_residence_time);
            }
          }
        }
      }

      if (min_residence_time != std::numeric_limits<double>::max() &&
          min_residence_time != 0.)
      {

        delta_time = min_residence_time / 100.;
      }
      else
      {
        // should throw  cause if delta <=0 we have incorrect timstep if first
        // branch fails delta_time is unchanged.
        // + If min_element fails flowmap might be invalid then default value is
        // not needed
        throw std::invalid_argument("No time step given and impossibe to "
                                    "estimate it with given flowmap");
      }
    }

    return delta_time;

    // auto min_liquid_residen_time = [](const auto& state) -> double
    // {
    //   constexpr double init_residence_time =
    //   std::numeric_limits<double>::max(); double min_residence_time =
    //   init_residence_time;

    //   for (size_t i = 0; i < state.n_compartments; ++i)
    //   {
    //     const auto view_flows = state.liquid_flow.getViewFlows();
    //     double sum_flows = 0.;
    //     for (size_t j = 0; j < state.n_compartments; ++j)
    //     {
    //       sum_flows += view_flows(i, j);
    //     }
    //     if (sum_flows != 0)
    //     {
    //       double residence_time = state.liquidVolume[i] / sum_flows;
    //       if (residence_time < min_residence_time && residence_time != 0.)
    //       {
    //         min_residence_time = residence_time;
    //       }
    //     }
    //   }

    //   if (min_residence_time <= 0 || min_residence_time ==
    //   init_residence_time)
    //   {
    //     throw std::invalid_argument("Flow map not valid");
    //   }

    //   return min_residence_time;
    // };

    // if (delta_time <= 0)
    // {
    //   auto minElement = std::min_element(
    //       iterator->begin(),
    //       iterator->end(),
    //       [&min_liquid_residen_time](const auto& state1, const auto& state2)
    //       {
    //         return min_liquid_residen_time(state1) <
    //                min_liquid_residen_time(state2);
    //       });

    //   if (minElement != iterator->end())
    //   {

    //     delta_time = min_liquid_residen_time(*minElement) / 100.;
    //   }
    //   else
    //   {
    //     // should throw  cause if delta <=0 we have incorrect timstep if
    //     first
    //     // branch fails delta_time is unchanged.
    //     // + If min_element fails flowmap might be invalid then default value
    //     is
    //     // not needed
    //     throw std::invalid_argument("No time step given and impossibe to "
    //                                 "estimate it with given flowmap");
    //   }
    // }

    // return delta_time;
  }

} // namespace

namespace Core
{
  GlobalInitialiser::GlobalInitialiser(const ExecInfo& _info,
                                       UserControlParameters _user_params,
                                       std::shared_ptr<IO::Logger> _logger)
      : info(_info), user_params(std::move(_user_params)),
        particle_per_process(0), is_host(info.current_rank == 0)

  {
    set_logger(std::move(_logger));

    Core::check_results_file_name(user_params);
    params = SimulationParameters::init(user_params);
    f_init_gas_flow =
        info.current_rank == 0 && params.is_two_phase_flow; // NOLINT
  }

  void GlobalInitialiser::set_logger(std::shared_ptr<IO::Logger> _logger)
  {
    if (_logger != nullptr)
    {
      this->logger = std::move(_logger);
    }
  }

  template <typename T> using OptionalPtr = GlobalInitialiser::OptionalPtr<T>;

  // CmtCommons::cma_exported_paths_t
  // GlobalInitialiser::get_path_files(const std::string& cma_case_path)
  // {
  //   CmtCommons::CMACaseInfo cma_case =
  //   CmtCommons::CMACaseInfoReader::load_case(
  //       cma_case_path + std::string(AutoGenerated::expected_cma_case_name));
  //   t_per_flowmap = cma_case.time_per_flowmap;
  //   return cma_case.paths;
  // }

  bool GlobalInitialiser::init_feed(
      std::optional<Simulation::Feed::SimulationFeed> _feed)
  {
    if (!check_steps(InitStep::InitState))
    {
      return false;
    }
    if (_feed) // TODO Improve error handling
    {
      auto index_max_compartments = liquid_volume.size() - 1;
      for (auto&& i : _feed->liquid_feeds())
      {
        // input_position: unsigned long so >0 no need to test this

        if (i.output_position > index_max_compartments)
        {
          return false;
        }
      }
      for (auto&& i : _feed->gas_feeds())
      {
        // input_position: unsigned long so >0 no need to test this

        if (i.output_position > index_max_compartments)
        {
          return false;
        }
      }
    }
    validate_step(InitStep::Feed);
    this->feed = std::move(_feed);
    return true;
  }

  std::optional<bool> GlobalInitialiser::init_state(
      const CmaUtils::TransitionnerPtrType& transiionner)
  {
    if (!check_steps(InitStep::Transitioner))
    {
      return std::nullopt;
    }
    if (is_host)
    {
      if (!host_init_state(transiionner).has_value())
      {
        return std::nullopt;
      }
    }
    mpi_broadcast();
    // liquid_neighbors.set_row_major();
    validate_step(InitStep::InitState);
    return true;
  }

  std::optional<CmaUtils::TransitionnerPtrType>
  GlobalInitialiser::init_transitionner()
  {

    auto d_transition = (is_host)
                            ? get_dtransitionner(user_params.cma_case_path)
                            : CmaUtils::TransitionnerPtrType::from_raw(nullptr);

    if (is_host)
    {
      // init_state(flow_handle);
      // Calculate the total number of time steps
      f_init_gas_flow = info.current_rank == 0 && params.is_two_phase_flow;
      const auto n_t =
          static_cast<size_t>(user_params.final_time / params.d_t) + 1;
      auto state = d_transition->get_current();
      if (logger)
      {
        // Note final "s"
        const auto str = (d_transition->size() > 1)
                             ? std::string(" flowmaps loaded with ")
                             : std::string(" flowmap loaded with ");
        const auto compartment_str = (state->n_compartments() > 1)
                                         ? std::string(" compartments")
                                         : std::string(" compartment");
        logger->print("Initializer",
                      IO::format(std::to_string(d_transition->size()),
                                 str,
                                 std::to_string(state->n_compartments()),
                                 compartment_str));
      }
    }

    validate_step(InitStep::Transitioner);
    init_state(d_transition);
    return d_transition;
  }

  OptionalPtr<Simulation::SimulationUnit> GlobalInitialiser::init_simulation(
      std::optional<Core::ScalarFactory::ScalarVariant> variant)
  {
    if (!check_steps(InitStep::Feed))
    {
      return std::nullopt;
    }
    auto scalar =
        variant.has_value() ? init_scalar(std::move(*variant)) : init_scalar();
    auto mc = init_monte_carlo();
    if (scalar.has_value() and mc.has_value())
    {

      return init_simulation(std::move(*mc), *scalar);
    }
    return std::nullopt;
  }

  OptionalPtr<MC::MonteCarloUnit> GlobalInitialiser::init_monte_carlo()
  {
    if (!check_steps(InitStep::InitState))
    {
      return std::nullopt;
    }

    const auto i_model =
        AutoGenerated::get_model_index_from_name(user_params.model_name);

    if (i_model == -2)
    {
      if (logger)
      {
        logger->alert("Initializer",
                      "Model not found, using DefaultModel instead");
      }
    }

    double total_mass = 0.;
    // TODO Add this as user param
    std::unique_ptr<ILoadBalancer> lb =
        std::make_unique<UniformLoadBalancer>(info.n_rank);

    const uint64_t particle_per_process =
        lb->balance(info.current_rank, user_params.number_particle);

    if (particle_per_process == 0)
    {
      throw std::runtime_error(
          "ERROR particle_per_process after balancing is 0");
    }

    if (particle_per_process > AutoGenerated::MC::MC_MAX_PARTICLE_BUFFER)
    {
      return std::nullopt;
    }

    auto mc_unit =
        AutoGenerated::wrap_init_model_selector(logger,
                                                i_model,
                                                particle_per_process,
                                                liquid_volume,
                                                flat_neighobrs,
                                                params.uniform_mc_init,
                                                total_mass);

    if (mc_unit == nullptr)
    {
      return std::nullopt;
    }

    double correct_total_mass = total_mass;
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
      correct_total_mass = 0;
#ifndef NO_MPI
      correct_total_mass = WrapMPI::all_reduce(total_mass);
#endif
    }

    MC::post_init_weight(
        mc_unit, user_params.biomass_initial_concentration, correct_total_mass);

    validate_step(InitStep::MC);

    return mc_unit;
  }

  OptionalPtr<Simulation::SimulationUnit> GlobalInitialiser::init_simulation(
      std::unique_ptr<MC::MonteCarloUnit> _unit,
      const Simulation::ScalarInitializer& scalar_init)
  {
    if (!check_steps(
            InitStep::Transitioner, InitStep::InitState, InitStep::Feed))
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

    // std::vector<double> kla(scalar_init.n_species);
    // kla[1] = 0.2; // 700 h-1
    // const auto type = Simulation::MassTransfer::Type::FixedKla{kla};

    auto simulation = std::make_unique<Simulation::SimulationUnit>(
        std::move(_unit), scalar_init, std::move(feed));

    // // FIXME
    // if (AutoGenerated::get_model_index_from_name(user_params.model_name) == 4
    // &&
    //     simulation->getDimensions().n_species < 4)
    // {
    //   throw std::runtime_error("Model must have 4 species concentrations");
    // }
    validate_step(InitStep::SimulationUnit);

    return simulation;
  }

  std::optional<bool> GlobalInitialiser::init_mtr_model(
      Simulation::SimulationUnit& unit,
      Simulation::MassTransfer::Type::MtrTypeVariant&& variant)
  {
    if (!check_steps(InitStep::SimulationUnit))
    {
      VERBOSE_ERROR
      return std::nullopt;
    }
    unit.setMtrModel(std::move(variant));

    // TODO
    return true;
  }

  std::optional<Simulation::ScalarInitializer>
  GlobalInitialiser::init_scalar(Core::ScalarFactory::ScalarVariant&& variant)
  {

    auto scalar_init = Core::ScalarFactory::scalar_factory(
        f_init_gas_flow, gas_volume, liquid_volume, variant);
    if (info.current_rank != 0) // FIXME
    {
      scalar_init.gas_flow = false;
    }
    validate_step(InitStep::Scalar);
    return scalar_init;
  }

  std::optional<Simulation::ScalarInitializer> GlobalInitialiser::init_scalar()
  {
    if (!check_steps(InitStep::Transitioner, InitStep::InitState))
    {
      return std::nullopt;
    }

    // We init gas flow specific object only on the host and if user defined two
    // phase flow simulation
    f_init_gas_flow = info.current_rank == 0 && params.is_two_phase_flow;

    Core::ScalarFactory::ScalarVariant arg;
    if (user_params.initialiser_path.empty())
    {
      if (logger)
      {
        logger->alert("Intializer", "Using Default scalar initialiser");
      }

      arg = params.is_two_phase_flow ? DefaultIntialiserTPF : DefaultIntialiser;
    }
    else
    {
      arg = Core::ScalarFactory::File(params.n_compartments,
                                      user_params.initialiser_path);
    }

    Simulation::ScalarInitializer scalar_init =
        Core::ScalarFactory::scalar_factory(
            f_init_gas_flow, gas_volume, liquid_volume, arg);

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
        if (logger)
        {
          logger->error(IO::format(
              "Initialisation step ", std::to_string(i), " not completed!"));
        }

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

  std::optional<bool> GlobalInitialiser::host_init_state(
      const CmaUtils::TransitionnerPtrType& transiionner)
  {
    // if (flow_handle == nullptr)
    // {
    //   throw std::runtime_error("Reactor not correclty initialised");
    // }
    // const auto* fstate = &flow_handle->get_unchecked(0);
    const auto& fstate = transiionner->get_current();

    params.n_compartments = fstate->n_compartments();
    // liquid_neighbors = fstate->f

    params.is_two_phase_flow = fstate->has_gas();

    auto slice_neighbors = fstate->flat_neighobrs();
    flat_neighobrs.resize(slice_neighbors.size());
    std::copy(
        slice_neighbors.begin(), slice_neighbors.end(), flat_neighobrs.begin());

    // worker_neighbor_data =
    // std::vector<size_t>(liquid_neighbors.data().begin(),
    //                                            liquid_neighbors.data().end());

    params.n_different_maps = transiionner->size();
    const auto lspan = fstate->get_liquid()->volume();
    if (lspan.empty())
    {
      return std::nullopt;
    }

    liquid_volume.clear();
    liquid_volume.resize(lspan.size());
    std::copy(lspan.begin(), lspan.end(), liquid_volume.begin());
    gas_volume.clear();
    gas_volume.resize(lspan.size());
    if (fstate->has_gas())
    {
      auto gspan = fstate->get_gas()->volume();

      std::copy(gspan.begin(), gspan.end(), gas_volume.begin());
    }

    if (liquid_volume.size() > 1)
    {
      params.d_t = get_time_step(user_params.delta_time, transiionner);
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
            : static_cast<size_t>(t_per_flowmap / params.d_t) + 1;

    params.n_per_flowmap = n_per_flowmap;

    return true;
  }

  void GlobalInitialiser::mpi_broadcast() // NOLINT
  {
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
#ifndef NO_MPI
      const std::size_t root_rank = 0;
      // Check once if we can broadcast integer, if not quit
      if (WrapMPI::broadcast(params.n_different_maps, root_rank) != 0)
      {
        WrapMPI::critical_error();
      }
      // broadcasting vector (malloc) is handled by WrapMPI so don't need to
      // resize before
      WrapMPI::broadcast(params.d_t, root_rank);
      WrapMPI::broadcast(params.n_per_flowmap, root_rank);
      WrapMPI::broadcast(params.n_compartments, root_rank);
      WrapMPI::broadcast(params.is_two_phase_flow, root_rank);
      WrapMPI::broadcast(liquid_volume, root_rank, info.current_rank);
      WrapMPI::broadcast(gas_volume, root_rank, info.current_rank);
      WrapMPI::broadcast(flat_neighobrs, root_rank, info.current_rank);
      // WrapMPI::bcst_iterator(_flow_handle, info.current_rank);
      if (!is_host)
      {
        // Neighbors data is a matrix with nrow the number of compartment (known
        // by host and workers) n_col is not directly known but the matrix size
        // is so we calculate column this way
        // assert(slice_neighbors.size() > params.n_compartments);
        // const size_t n_col = slice_neighbors.size() / params.n_compartments;
        // assert(slice_neighbors.size() % params.n_compartments == 0);

        // liquid_neighbors = CmaRead::Neighbors::Neighbors_const_view_t(
        //     slice_neighbors, params.n_compartments, n_col);
      }
      if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
      {
        WrapMPI::barrier(); // This barrier is probably useless
      }

#endif
    }
  }

} // namespace Core
