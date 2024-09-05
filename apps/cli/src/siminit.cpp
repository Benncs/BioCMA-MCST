#include "scalar_factory.hpp"
#include <cassert>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/neighbors.hpp>
#include <cmt_common/cma_case.hpp>
#include <common/common.hpp>
#include <common/simulation_parameters.hpp>
#include <core/cp_flag.hpp>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <mc/mcinit.hpp>
#include <mc/prng/distribution.hpp>
#include <memory>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <wrap_init_model_selector.hpp>

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#else
namespace MPI_W
{
  void barrier();
};
#endif

static const ScalarFactory::ScalarVariant DefaultIntialiser =
    ScalarFactory::Uniform({1.});

/**
 * @brief Initializes and retrieves the reactor state based on provided
 * parameters and flow handle.
 *
 * This function read and returns a pointer to a constant
 * `CmaRead::ReactorState` object. The state is read using the given
 * simulation parameters, flow iterator, and exported paths.
 *
 * @param params The `SimulationParameters` object containing settings for
 * initialization.
 * @param flow_handle A unique pointer to the `CmaRead::FlowIterator` for
 * to read flowmap.
 * @param paths The `CmtCommons::cma_exported_paths_t` containing paths needed
 * for initialization.
 * @return A pointer to a constant `CmaRead::ReactorState` representing the
 * initialized state.
 */
static CmaRead::ReactorState const *
init_state(SimulationParameters &params,
           std::unique_ptr<CmaRead::FlowIterator> &flow_handle,
           const CmtCommons::cma_exported_paths_t &paths);

/**
 * @brief Initialisation step required only for the host rank
 */
static void
init_host_only(const ExecInfo &info,
               SimulationParameters &params,
               std::unique_ptr<CmaRead::FlowIterator> &_flow_handle,
               CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors,
               std::vector<double> &liq_volume,
               std::vector<double> &gas_volume,
               std::vector<size_t> &worker_neighbor_data);

std::unique_ptr<Simulation::SimulationUnit>
init_simulation(const ExecInfo &info,
                SimulationParameters &params,
                std::unique_ptr<Simulation::FlowMapTransitioner> &transitioner)
{

  const auto &user_params = params.user_params;

  // Declaration of variable that will be send between Host and Worker
  // In case of non multirank context some of the following operation are
  // useless but they are idempotent so this will not interfere with data
  // integrity
  std::vector<double> liquid_volume;
  std::vector<double> gas_volume;
  std::vector<size_t> worker_neighbor_data;
  CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors;
  std::unique_ptr<CmaRead::FlowIterator> _flow_handle = nullptr;

  // Host is namely going to read flowmap from local files,
  // From ReactorState is read from files and the previous data is filled
  init_host_only(info,
                 params,
                 _flow_handle,
                 liquid_neighbors,
                 liquid_volume,
                 gas_volume,
                 worker_neighbor_data);

  // After data being read, host send reactor and simulation data to workers
  if constexpr (FlagCompileTIme::use_mpi)
  {
#ifndef NO_MPI
    // Check once if we can broadcast integer, if not quit
    if (MPI_W::broadcast(params.n_different_maps, 0) != 0)
    {
      MPI_W::critical_error();
    }
    // broadcasting vector (malloc) is handled by MPI_W so don't need to resize
    // before
    MPI_W::broadcast(params.d_t, 0);
    MPI_W::broadcast(params.n_per_flowmap, 0);
    MPI_W::broadcast(params.n_compartments, 0);
    MPI_W::broadcast(params.is_two_phase_flow, 0);
    MPI_W::broadcast(liquid_volume, 0, info.current_rank);
    MPI_W::broadcast(gas_volume, 0, info.current_rank);
    MPI_W::broadcast(worker_neighbor_data, 0, info.current_rank);
#endif
  }

  // MPI_W::bcst_iterator(_flow_handle, info.current_rank);
  if (info.current_rank != 0)
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
  liquid_neighbors.set_row_major();

  // MC::UniformLawINT law_param = {0,
  //                                static_cast<int>(params.n_compartments -
  //                                1)};

  if constexpr (FlagCompileTIme::use_mpi)
  {
    MPI_W::barrier(); // This barrier is probably useless
  }

  // MonteCarlo Initialisation
  auto mc_unit = AutoGenerated::wrap_init_model_selector(
      AutoGenerated::get_model_index_from_name(user_params.model_name),
      info,
      user_params.numper_particle,
      liquid_volume,
      liquid_neighbors,
      user_params.biomass_initial_concentration);

  // We init gas flow specific object only on the host and if user defined two
  // phase flow simulation
  const bool f_init_gas_flow =
      info.current_rank == 0 && params.is_two_phase_flow;

  // std::vector<double> concentrations = {1, 5};
  // std::vector<size_t> indices = {0, 5, 1500};
  // ScalarFactory::Local arg = {concentrations, indices};
  // ScalarFactory::Uniform uniform = {concentrations, std::nullopt};
  ScalarFactory::ScalarVariant arg;
  if (params.user_params.initialiser_path == "")
  {
    std::cerr << "WARNING: using Default Initialiser" << std::endl;
    arg = DefaultIntialiser;
  }
  else
  {
    arg = ScalarFactory::File(params.n_compartments,
                              params.user_params.initialiser_path);
  }

  auto scalar_init = ScalarFactory::scalar_factory(
      f_init_gas_flow, gas_volume, liquid_volume, arg);

  // Construct the main simulation object (one per rank)
  auto simulation = std::make_unique<Simulation::SimulationUnit>(
      std::move(mc_unit), scalar_init);

  // Calculate the total number of time steps
  const auto n_t = static_cast<size_t>(user_params.final_time / params.d_t) + 1;

  // Transitioner handles flowmap transition between time step, flowmaps are
  // only located in host but transitioner handles cache and receiving for
  // workers
  transitioner = std::make_unique<Simulation::FlowMapTransitioner>(
      params.n_different_maps,
      params.n_per_flowmap,
      Simulation::FlowMapTransitioner::Discontinuous,
      n_t,
      std::move(_flow_handle),
      f_init_gas_flow);

  std::cout << "INIT OK" << std::endl;
  return simulation;
}

static void
init_host_only(const ExecInfo &info,
               SimulationParameters &params,
               std::unique_ptr<CmaRead::FlowIterator> &_flow_handle,
               CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors,
               std::vector<double> &liq_volume,
               std::vector<double> &gas_volume,
               std::vector<size_t> &worker_neighbor_data)
{

  if (info.current_rank != 0)
  {
    return;
  }

  // Read CMA data
  std::string case_name = params.user_params.cma_case_path + "/cma_case";
  CmtCommons::CMACaseInfo cma_case =
      CmtCommons::CMACaseInfoReader::load_case(case_name);

  // Init _flow_handle and return the first state needed to get some flow
  // information All states have the same properties (number of compartment,
  // number of neighbors ...)
  const CmaRead::ReactorState *fstate =
      init_state(params, _flow_handle, cma_case.paths);

  params.n_compartments = fstate->n_compartments;
  liquid_neighbors = fstate->liquid_flow.getViewNeighors();
  if (fstate->gas_flow.is_empty())
  {
    params.is_two_phase_flow = false;
  }

  worker_neighbor_data = std::vector<size_t>(liquid_neighbors.data().begin(),
                                             liquid_neighbors.data().end());

  params.n_different_maps = _flow_handle->size();

  liq_volume = fstate->liquidVolume;
  gas_volume = fstate->gasVolume;

  if (liq_volume.size() > 1)
  {
    // When multiple compartments are present, it implies the existence of
    // internal hydrodynamic time scales. To account for this, the simulation's
    // explicit time step is calculated to approximate a CFL condition, with the
    // formula: time_step = min(residence_time) / 100. This approach ensures
    // that the fluid movement between two steps is accurately represented
    // without losing flow information.

    params.d_t = (params.user_params.delta_time == 0.)
                     ? _flow_handle->MinLiquidResidenceTime() / 100.
                     : params.user_params.delta_time;
  }
  else
  {
    // IF 0D reactor, no hydrodynamic time scale (internal)
    params.d_t = params.user_params.delta_time;
  }

  // Define the duration of each flowmap and compute steps per flowmap
  const double t_per_flowmap = cma_case.time_per_flowmap;
  params.t_per_flow_map = t_per_flowmap;
  const auto n_per_flowmap =
      (t_per_flowmap == 0 || params.n_different_maps == 1)
          ? 1
          : static_cast<size_t>(t_per_flowmap /
                                static_cast<double>(params.d_t)) +
                1;

  params.n_per_flowmap = n_per_flowmap;

  register_run(info, params);
}

static CmaRead::ReactorState const *
init_state(SimulationParameters &params,
           std::unique_ptr<CmaRead::FlowIterator> &flow_handle,
           const CmtCommons::cma_exported_paths_t &paths)
{
  // Ensure that state is correcly loaded
  CmaRead::ReactorState const *state = nullptr;
  try
  {
    flow_handle =
        std::make_unique<CmaRead::FlowIterator>(params.flow_files, paths);
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
    std::cout << "Flowmap loaded with " << state->n_compartments
              << " compartments" << std::endl;
  }
  catch (const std::exception &e)
  {
    std::stringstream err;
    err << "Error while reading files\t:" << e.what();
    throw std::runtime_error(err.str());
  }

  return state;
}
