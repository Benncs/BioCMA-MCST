#include "cma_read/light_2d_view.hpp"
#include "cma_read/neighbors.hpp"
#include "models/types.hpp"
#include "mpi.h"
#include "mpi_w/impl_op.hpp"
#include "mpi_w/message_t.hpp"
#include <cstddef>
#include <memory>
#include <siminit.hpp>

#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <cstdio>
#include <fstream>
#include <mc/mcinit.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <rt_init.hpp>
#include <simulation/simulation.hpp>
#include <stdexcept>
#include <utility>
#include <vector>
static ReactorState const *
init_state(SimulationParameters &params,
           std::shared_ptr<FlowIterator> &flow_handle);
static void register_run(ExecInfo &exec, SimulationParameters &params);

Simulation::SimulationUnit
init_simulation(ExecInfo &info,
                SimulationParameters &params,
                std::shared_ptr<FlowIterator> &_flow_handle,
                KModel model)
{

  std::vector<double> liq_volume;
  std::vector<double> gas_volume;
  std::vector<size_t> worker_neighbor_data;

  size_t *worker_neighbor_data_ptr = nullptr;

  size_t worker_neighbor_data_size = 0;

  Neighbors::Neighbors_const_view_t liquid_neighbors;
  size_t n_compartments = 0;
  if (info.current_rank == 0)
  {
    const ReactorState *fstate = init_state(params, _flow_handle);

    n_compartments = fstate->n_compartments;
    liquid_neighbors = fstate->liquid_flow.getViewNeighors();

    worker_neighbor_data_size = liquid_neighbors.size();
    worker_neighbor_data_ptr = const_cast<size_t *>(
        liquid_neighbors.data()
            .data()); // TODO BIREM add rawdata() to get ptr directly

    params.n_different_maps = _flow_handle->loop_size();
    liq_volume = fstate->liquidVolume;
    gas_volume = fstate->gasVolume;
    params.d_t = (params.d_t == 0.)
                     ? _flow_handle->MinLiquidResidenceTime() / 100.
                     : params.d_t;
    auto n_t = static_cast<size_t>(params.final_time / params.d_t);
    _flow_handle->setRepetition(n_t / params.n_different_maps);
    register_run(info, params);
  }

  if (MPI_W::broadcast(params.n_different_maps, 0) != 0)
  {
    MPI_W::critical_error();
  }

  MPI_W::broadcast(params.d_t, 0);
  MPI_W::broadcast(worker_neighbor_data_size, 0);
  MPI_W::broadcast(liq_volume, 0, info.current_rank);
  MPI_W::broadcast(gas_volume, 0, info.current_rank);

  if (info.current_rank != 0)
  {
    n_compartments = liq_volume.size();
    worker_neighbor_data.resize(worker_neighbor_data_size, 0);
    worker_neighbor_data_ptr = worker_neighbor_data.data();
  }

  if (MPI_W::_broadcast_unsafe(worker_neighbor_data_ptr,
                               worker_neighbor_data_size,
                               0) != MPI_SUCCESS)
  {
    MPI_W::critical_error();
  }

  if (info.current_rank != 0)
  {
    auto n_col = worker_neighbor_data_size / n_compartments;
    liquid_neighbors = L2DView<const size_t>(
        worker_neighbor_data, n_compartments, n_col, false);
  }

  MPI_W::barrier();
  auto unit = MC::init_unit(info, liq_volume, liquid_neighbors);

  auto container = MC::init_container(info, params.n_particles);

  // size_t unit_size =
  // unit->domain.getNumberCompartments()*sizeof(unit->domain[0]) +
  // sizeof(unit->domain.getNeighbors(0))*unit->domain.getNeighbors(0).size()*unit->domain.getNeighbors().size();
  // size_t container_size =
  // container->to_process.size()*sizeof(container->to_process[0]);

  auto simulation = Simulation::SimulationUnit(info,
                                               std::move(unit),
                                               std::move(container),
                                               params.n_species,
                                               model,
                                               info.current_rank == 0);
  return simulation;
}

static ReactorState const *
init_state(SimulationParameters &params,
           std::shared_ptr<FlowIterator> &flow_handle)
{
  ReactorState const *state = nullptr;
  try
  {
    flow_handle = std::make_shared<FlowIterator>(params.flow_files);
    if (flow_handle == nullptr)
    {
      throw std::runtime_error("Flow map are not loaded");
    }
    std::cout << "Flowmap loaded: " << flow_handle->loop_size() << std::endl;

    state = &flow_handle->operator()(0);
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

static void register_run(ExecInfo &exec, SimulationParameters &params)
{
  // Open the file in append mode
  std::ofstream env(env_file_path(), std::ios_base::app);
  if (env.is_open())
  {
    append_date_time(env);
    env << exec;
    env << params;
    env << std::endl;
  }
  else
  {
    std::cerr << "Error: Unable to open file for writing\n";
  }
}