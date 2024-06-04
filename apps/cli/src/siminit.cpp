

#include "mpi_w/iteration_payload.hpp"
#include <cstddef>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <birem_common/cma_case.hpp>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/neighbors.hpp>
#include <common/common.hpp>
#include <mc/mcinit.hpp>
#include <mc/prng/distribution.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <simulation/simulation.hpp>

#include <simulation/update_flows.hpp>

static ReactorState const *
init_state(SimulationParameters &params,
           std::unique_ptr<FlowIterator> &flow_handle,
           const cma_exported_paths_t &paths);

static void init_host_only(ExecInfo &info,
                           SimulationParameters &params,
                           std::unique_ptr<FlowIterator> &_flow_handle,
                           Neighbors::Neighbors_const_view_t liquid_neighbors,
                           std::vector<double> &liq_volume,
                           std::vector<double> &gas_volume,
                           std::vector<size_t> &worker_neighbor_data);

Simulation::SimulationUnit
init_simulation(ExecInfo &info,
                SimulationParameters &params,

                std::unique_ptr<Simulation::FlowMapTransitioner> &transitioner,
                KModel model,
                MC::DistributionVariantInt &&initial_particle_distribution)
{

  std::vector<double> liq_volume;
  std::vector<double> gas_volume;
  std::vector<size_t> worker_neighbor_data;
  Neighbors::Neighbors_const_view_t liquid_neighbors;
  std::unique_ptr<FlowIterator> _flow_handle = nullptr;
  init_host_only(info,
                 params,
                 _flow_handle,
                 liquid_neighbors,
                 liq_volume,
                 gas_volume,
                 worker_neighbor_data);

  if (MPI_W::broadcast(params.n_different_maps, 0) != 0)
  {
    MPI_W::critical_error();
  }

  MPI_W::broadcast(params.d_t, 0);
  MPI_W::broadcast(params.n_per_flowmap, 0);
  MPI_W::broadcast(params.n_compartments, 0);
  MPI_W::broadcast(params.is_two_phase_flow, 0);
  MPI_W::broadcast(liq_volume, 0, info.current_rank);
  MPI_W::broadcast(gas_volume, 0, info.current_rank);
  MPI_W::broadcast(worker_neighbor_data, 0, info.current_rank);

  // MPI_W::bcst_iterator(_flow_handle, info.current_rank);
  if (info.current_rank != 0)
  {

    auto n_col = worker_neighbor_data.size() / params.n_compartments;
    liquid_neighbors = L2DView<const size_t>(
        worker_neighbor_data, params.n_compartments, n_col, false);
  }

  MPI_W::barrier();

  auto mc_unit =
      MC::init(info, params.n_particles, liq_volume, liquid_neighbors);

  bool tpf = info.current_rank == 0 && params.is_two_phase_flow;

  auto simulation =
      Simulation::SimulationUnit(info,
                                 std::move(mc_unit),
                                 gas_volume,
                                 liq_volume,
                                 params.n_species,
                                 model,
                                 std::move(initial_particle_distribution),
                                 tpf);
  // Calculate the total number of time steps
  const auto n_t = static_cast<size_t>(params.final_time / params.d_t) + 1;

  transitioner = std::make_unique<Simulation::FlowMapTransitioner>(
      params.n_different_maps,
      params.n_per_flowmap,
      Simulation::FlowMapTransitioner::Discontinuous,
      n_t,
      std::move(_flow_handle),
      params.is_two_phase_flow);

  return simulation;
}

static void init_host_only(ExecInfo &info,
                           SimulationParameters &params,
                           std::unique_ptr<FlowIterator> &_flow_handle,
                           Neighbors::Neighbors_const_view_t liquid_neighbors,
                           std::vector<double> &liq_volume,
                           std::vector<double> &gas_volume,
                           std::vector<size_t> &worker_neighbor_data)
{

  if (info.current_rank != 0)
  {
    return;
  }
  CMACaseInfo cma_case =
      CMACaseInfoReader::load_case(params.root + "/cma_case");

  const ReactorState *fstate = init_state(params, _flow_handle, cma_case.paths);

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
    params.d_t = (params.d_t == 0.)
                     ? _flow_handle->MinLiquidResidenceTime() / 100.
                     : params.d_t;
  }

  // const auto n_t = static_cast<size_t>(params.final_time / params.d_t) + 1;

  // Define the duration of each flowmap and compute steps per flowmap
  const double t_per_flowmap = 0;//cma_case.time_per_flowmap;

  const auto n_per_flowmap =
      (t_per_flowmap == 0 || params.n_different_maps == 1)
          ? 1
          : static_cast<size_t>(t_per_flowmap /
                                static_cast<double>(params.d_t)) +
                1;

  // // Calculate the number of repetitions: Not needed anymore
  // const size_t n_repetition =
  //     n_t / (params.n_different_maps * n_per_flowmap) + 1;

  params.n_per_flowmap = n_per_flowmap;

  // _flow_handle->setRepetition(n_repetition, n_per_flowmap);
  register_run(info, params);
}

static ReactorState const *
init_state(SimulationParameters &params,
           std::unique_ptr<FlowIterator> &flow_handle,
           const cma_exported_paths_t &paths)
{
  ReactorState const *state = nullptr;
  try
  {
    flow_handle = std::make_unique<FlowIterator>(params.flow_files, paths);
    if (flow_handle == nullptr)
    {
      throw std::runtime_error("Flow map are not loaded");
    }
    std::cout << "Flowmap loaded: " << flow_handle->size() << std::endl;

    state = &flow_handle->get_unchecked(0);
  }
  catch (const std::exception &e)
  {
    throw std::runtime_error("Error while reading files");
  }

  if (state == nullptr)
  {
    throw std::runtime_error("Error while reading files");
  }
  return state;
}
