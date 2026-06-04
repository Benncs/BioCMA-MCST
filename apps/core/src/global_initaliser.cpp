#include <algorithm>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <cma_utils/alias.hpp>
#include <common/env_var.hpp>
#include <common/logger.hpp>
#include <core/global_initaliser.hpp>
#include <core/scalar_factory.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <cstdio>
#include <dataexporter/data_exporter.hpp>
#include <exception>
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
  const Core::ScalarFactory::ScalarVariant DefaultIntialiser
      = Core::ScalarFactory::Uniform({ 1. });

  const Core::ScalarFactory::ScalarVariant DefaultIntialiserTPF
      = Core::ScalarFactory::Local({ 0., 0. },
                                   { 0 },
                                   std::vector<double>({ 0, 0. }),
                                   std::vector<size_t>({ 0 }));

  double
  get_time_step(double user_deta_time,
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

      const auto min_residence_time
          = CmaUtils::get_min_residence_time(iterator);

      if (min_residence_time != std::numeric_limits<double>::max()
          && min_residence_time != 0.)
      {

        delta_time = min_residence_time / 100.;
      }
      else
      {
        // should throw  cause if delta <=0 we have incorrect timstep if
        // first branch fails delta_time is unchanged.
        // + If min_element fails flowmap might be invalid then default value
        // is not needed
        throw std::invalid_argument("No time step given and impossibe to "
                                    "estimate it with given flowmap");
      }
    }

    return delta_time;
  }

  // FIXME
  std::unique_ptr<ILoadBalancer>
  lb_factory(uint32_t s)
  {
    auto bounded = Common::read_env<uint32_t>("BIOMC_LBBOUND");
    if (bounded)
    {
      return std::make_unique<BoundLoadBalancer>(s, *bounded);
    }
    else
    {
      return std::make_unique<UniformLoadBalancer>(s);
    }
  }

  size_t
  compute_n_per_flowmap(double t_per_flowmap,
                        size_t n_different_maps,
                        double d_t)
  {
    if (t_per_flowmap == 0 || n_different_maps == 1)
    {
      return 1;
    }
    // TODO: fix float/integer division
    return static_cast<size_t>(t_per_flowmap / d_t) + 1;
  }

} // namespace

namespace Core
{
  template <typename T> using OptionalPtr = GlobalInitialiser::OptionalPtr<T>;
  GlobalInitialiser::GlobalInitialiser(const ExecInfo& _info,
                                       UserControlParameters _user_params,
                                       std::shared_ptr<IO::Logger> _logger)
      : info(_info), user_params(std::move(_user_params)),
        is_host(info.current_rank == 0)

  {
    set_logger(std::move(_logger));

    if (!Core::fill_and_check_result_file_path(this->m_logger, user_params))
    {
      throw std::runtime_error("TODO bad path");
    };
    m_params = SimulationParameters::init(user_params);
    f_init_gas_flow
        = info.current_rank == 0 && m_params.is_two_phase_flow; // NOLINT
  }

  void
  GlobalInitialiser::set_logger(std::shared_ptr<IO::Logger> _logger)
  {
    if (_logger != nullptr)
    {
      this->m_logger = std::move(_logger);
    }
  }

  bool
  GlobalInitialiser::init_feed(
      std::optional<Simulation::Feed::SimulationFeed> _feed)
  {

    if (!check_steps(InitStep::InitState))
    {
      return false;
    }
    if (_feed) // TODO Improve error handling
    {
      auto index_max_compartments = m_liquid_volume.size() - 1;

      auto pred = [&index_max_compartments](const auto& i)
      { return i.output_position <= index_max_compartments; };

      if (!std::ranges::all_of(_feed->liquid_feeds(), pred))
      {
        return false;
      }

      if (!std::ranges::all_of(_feed->gas_feeds(), pred))
      {
        return false;
      }
    }
    validate_step(InitStep::Feed);
    this->m_feed = std::move(_feed);
    return true;
  }

  std::optional<bool>
  GlobalInitialiser::init_state(
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
    auto d_transition = CmaUtils::TransitionnerPtrType::from_raw(nullptr);

    if (is_host)
    {
      try
      {
        d_transition = (is_host)
                           ? get_dtransitioner(user_params.cma_case_path)
                           : CmaUtils::TransitionnerPtrType::from_raw(nullptr);
      }
      catch (std::exception& e)
      {
        if (m_logger)
        {
          m_logger->error(e.what());
        }
        return std::nullopt;
      }

      const auto state = d_transition->get_current();
      f_init_gas_flow = info.current_rank == 0 && state->has_gas();
      if (m_logger)
      {
        // Note final "s"
        const auto str = (d_transition->size() > 1)
                             ? std::string(" flowmaps loaded with ")
                             : std::string(" flowmap loaded with ");
        const auto compartment_str = (state->n_compartments() > 1)
                                         ? std::string(" compartments")
                                         : std::string(" compartment");
        m_logger->print("Initializer",
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

  OptionalPtr<MC::MonteCarloUnit>
  GlobalInitialiser::init_monte_carlo()
  {
    if (!check_steps(InitStep::InitState, InitStep::Feed))
    {
      return std::nullopt;
    }

    std::size_t n_samples = 0;
    if (this->m_feed)
    {
      n_samples += (this->m_feed->n_liquid_flow() != 0) ? 1 : 0;
    }
    n_samples += (this->m_liquid_volume.size() > 1) ? 2 : 0;

    const auto i_model
        = AutoGenerated::get_model_index_from_name(user_params.model_name);

    if (i_model == -2)
    {
      if (m_logger)
      {
        m_logger->alert("Initializer",
                        "Model not found, using DefaultModel instead");
      }
    }

    double total_mass = 0.;
    // TODO Add this as user param
    std::unique_ptr<ILoadBalancer> lb = lb_factory(info.n_rank);

    const uint64_t particle_per_process
        = lb->balance(info.current_rank, user_params.number_particle);

    if (particle_per_process == 0)
    {
      throw std::runtime_error(
          "ERROR particle_per_process after balancing is 0");
    }

    if (particle_per_process > AutoGenerated::MC::MC_MAX_PARTICLE_BUFFER)
    {
      return std::nullopt;
    }

    auto mc_unit
        = AutoGenerated::wrap_init_model_selector(m_logger,
                                                  i_model,
                                                  particle_per_process,
                                                  n_samples,
                                                  m_liquid_volume,
                                                  m_flat_neighobrs,
                                                  m_params.uniform_mc_init,
                                                  total_mass);

    if (mc_unit == nullptr)
    {
      return std::nullopt;
    }

    double correct_total_mass = total_mass;
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
      // correct_total_mass = 0; //
#ifndef NO_MPI
      correct_total_mass = WrapMPI::all_reduce(total_mass);
#else
      __builtin_unreachable(); // FIXME
      throw std::runtime_error("ERROR MPI activate but not implemented ");
#endif
    }

    if (correct_total_mass <= 0.)
    {
      throw std::runtime_error("ERROR model not found returns negative mass");
    }

    MC::post_init_weight(
        mc_unit, user_params.biomass_initial_concentration, correct_total_mass);

    validate_step(InitStep::MC);

    return mc_unit;
  }

  OptionalPtr<Simulation::SimulationUnit>
  GlobalInitialiser::init_simulation(
      std::optional<Core::ScalarFactory::ScalarVariant> variant)
  {
    if (!check_steps(InitStep::Feed))
    {
      return std::nullopt;
    }
    auto scalar = variant.has_value() ? init_scalar(std::move(*variant))
                                      : init_scalar();

    auto mc = init_monte_carlo();
    if (scalar.has_value() && mc.has_value())
    {

      if ((*mc)->getSpeciesNames().size() > (*scalar).n_species)
      {
        // throw std::invalid_argument("Should declare at least number of
        // species "
        //                             "required by chosen model");
        return std::nullopt;
      }

      // Start long forwarding of scalar_init over all subsequent calls have
      // && type and need to std::move(scalar)
      return init_simulation(std::move(*mc), std::move(*scalar));
    }
    return std::nullopt;
  }

  OptionalPtr<Simulation::SimulationUnit>
  GlobalInitialiser::init_simulation(
      std::unique_ptr<MC::MonteCarloUnit> _unit,
      Simulation::ScalarInitializer&& scalar_init)
  {
    if (!check_steps(InitStep::Transitioner,
                     InitStep::InitState,
                     InitStep::Feed,
                     InitStep::MC))
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

    // std::vector<double> kla(scalar_init.n_species);
    // kla[1] = 0.2; // 700 h-1
    // const auto type = Simulation::MassTransfer::Type::FixedKla{kla};

    auto simulation = std::make_unique<Simulation::SimulationUnit>(
        std::move(_unit), std::move(scalar_init), std::move(m_feed));

    // // FIXME
    // if (AutoGenerated::get_model_index_from_name(user_params.model_name) ==
    // 4
    // &&
    //     simulation->getDimensions().n_species < 4)
    // {
    //   throw std::runtime_error("Model must have 4 species concentrations");
    // }
    validate_step(InitStep::SimulationUnit);

    return simulation;
  }

  std::optional<bool>
  GlobalInitialiser::init_mtr_model(
      Simulation::SimulationUnit& unit,
      std::optional<Simulation::MassTransfer::Type::MtrTypeVariant>&& variant)
  {
    if (!check_steps(InitStep::SimulationUnit))
    {
      VERBOSE_ERROR
      return std::nullopt;
    }
    if (variant)
    {
      unit.setMtrModel(std::move(*variant));
    }

    validate_step(InitStep::MTR);
    // TODO
    return true;
  }

  std::optional<Simulation::ScalarInitializer>
  GlobalInitialiser::init_scalar(Core::ScalarFactory::ScalarVariant&& variant)
  {

    auto scalar_init = Core::ScalarFactory::scalar_factory(
        f_init_gas_flow, m_gas_volume, m_liquid_volume, variant);
    if (info.current_rank != 0) // FIXME
    {
      scalar_init.gas_flow = false;
    }
    validate_step(InitStep::Scalar);
    return scalar_init;
  }

  std::optional<Simulation::ScalarInitializer>
  GlobalInitialiser::init_scalar()
  {
    if (!check_steps(InitStep::Transitioner, InitStep::InitState))
    {
      return std::nullopt;
    }

    // We init gas flow specific object only on the host and if user defined
    // two phase flow simulation
    f_init_gas_flow = info.current_rank == 0 && m_params.is_two_phase_flow;

    Core::ScalarFactory::ScalarVariant arg;
    if (user_params.initialiser_path.empty())
    {
      if (m_logger)
      {
        m_logger->alert("Intializer", "Using Default scalar initialiser");
      }

      arg = m_params.is_two_phase_flow ? DefaultIntialiserTPF
                                       : DefaultIntialiser;
    }
    else
    {
      arg = Core::ScalarFactory::File(m_params.n_compartments,
                                      user_params.initialiser_path);
    }

    Simulation::ScalarInitializer scalar_init
        = Core::ScalarFactory::scalar_factory(
            f_init_gas_flow, m_gas_volume, m_liquid_volume, arg);

    if (info.current_rank != 0) // FIXME
    {
      scalar_init.gas_flow = false;
    }

    validate_step(InitStep::Scalar);

    return scalar_init;
  }

  void
  GlobalInitialiser::set_initial_number_particle(uint64_t np) noexcept
  {
    m_params.number_particle = np;
  }

  bool
  GlobalInitialiser::check_init_terminate() const
  {
    int i = 0;

    for (auto&& step : validated_steps)
    {
      if (!step)
      {
        if (m_logger)
        {
          m_logger->error(IO::format(
              "Initialisation step ", std::to_string(i), " not completed!"));
        }

        return false;
      }
      ++i;
    }

    return true;
  }

  [[nodiscard]] SimulationParameters
  GlobalInitialiser::get_parameters() const
  {

    return m_params;
  }

  std::optional<bool>
  GlobalInitialiser::host_init_state(
      const CmaUtils::TransitionnerPtrType& transiionner)
  {

    const auto& fstate = transiionner->get_current();
    m_params.n_compartments = fstate->n_compartments();
    m_params.is_two_phase_flow = fstate->has_gas();
    m_params.n_different_maps = transiionner->size();

    const auto& slice_neighbors = fstate->flat_neighobrs();
    m_flat_neighobrs.assign(slice_neighbors.begin(), slice_neighbors.end());

    // Volumes
    const auto liquid_volume_span = fstate->get_liquid()->volume();
    if (liquid_volume_span.empty())
    {
      return std::nullopt;
    }

    m_liquid_volume.assign(liquid_volume_span.begin(),
                           liquid_volume_span.end());
    m_gas_volume.assign(liquid_volume_span.size(), 0.0);

    if (fstate->has_gas())
    {
      const auto gas_volume_span = fstate->get_gas()->volume();
      m_gas_volume.assign(gas_volume_span.begin(), gas_volume_span.end());
    }

    if (m_liquid_volume.size() > 1)
    {
      m_params.d_t = get_time_step(user_params.delta_time, transiionner);
    }
    else
    {
      // IF 0D reactor, no hydrodynamic time scale (internal)
      m_params.d_t = user_params.delta_time;
    }

    // Define the duration of each flowmap and compute steps per flowmap
    m_params.t_per_flow_map = t_per_flowmap;

    m_params.n_per_flowmap = compute_n_per_flowmap(
        t_per_flowmap, m_params.n_different_maps, m_params.d_t);

    return true;
  }

  void
  GlobalInitialiser::mpi_broadcast() // NOLINT
  {
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
#ifndef NO_MPI
      const std::size_t root_rank = 0;
      // Check once if we can broadcast integer, if not quit
      if (WrapMPI::broadcast(m_params.n_different_maps, root_rank) != 0)
      {
        WrapMPI::critical_error();
      }
      // broadcasting vector (malloc) is handled by WrapMPI so don't need to
      // resize before
      WrapMPI::broadcast(m_params.d_t, root_rank);
      WrapMPI::broadcast(m_params.n_per_flowmap, root_rank);
      WrapMPI::broadcast(m_params.n_compartments, root_rank);
      WrapMPI::broadcast(m_params.is_two_phase_flow, root_rank);
      WrapMPI::broadcast(m_liquid_volume, root_rank, info.current_rank);
      WrapMPI::broadcast(m_gas_volume, root_rank, info.current_rank);
      WrapMPI::broadcast(m_flat_neighobrs, root_rank, info.current_rank);

      if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
      {
        WrapMPI::barrier(); // This barrier is probably useless
      }

#endif
    }
  }

} // namespace Core
