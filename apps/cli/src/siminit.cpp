#include "models/types.hpp"
#include <memory>
#include <siminit.hpp>

#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <cstdio>
#include <mpi_w/wrap_mpi.hpp>
#include <simulation/simulation.hpp>
#include <rt_init.hpp>
#include <mc/mcinit.hpp>
#include <stdexcept>
#include <utility>
#include <vector>
#include <fstream>
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
  std::vector<std::vector<size_t>> liquid_neighbors;
  size_t n_compartments = 0;
  size_t nmap = 0;
  double opti_dt = params.d_t;
  if (info.current_rank == 0)
  {
    const ReactorState *fstate = init_state(params, _flow_handle);

    n_compartments = fstate->n_compartments;
    liquid_neighbors = fstate->liquid_flow.neigbors;
    nmap = _flow_handle->loop_size();
    liq_volume = fstate->liquidVolume;
    gas_volume = fstate->gasVolume;
    opti_dt = (params.d_t == 0.) ? _flow_handle->MinLiquidResidenceTime() / 100.
                                 : params.d_t;
    auto n_t = static_cast<size_t>(params.final_time / opti_dt);
    _flow_handle->setRepetition(n_t / nmap);
    register_run(info, params);
  }

  // MPI_W::broadcast(d_t, 0);
  if (MPI_W::broadcast(nmap, 0) != 0)
  {
    MPI_W::critical_error();
  }

  MPI_W::broadcast(opti_dt, 0);
  params.d_t = opti_dt;

  MPI_W::broadcast(liq_volume, 0, info.current_rank);
  MPI_W::broadcast(gas_volume, 0, info.current_rank);
  params.n_different_maps = nmap;
  if (info.current_rank != 0)
  {
    n_compartments = liq_volume.size();
    liquid_neighbors.resize(n_compartments);
  }

  for (auto &&i : liquid_neighbors)
  {
    MPI_W::broadcast(i, 0, info.current_rank);
  }

  auto unit = MC::init_unit(info, liq_volume, std::move(liquid_neighbors));

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
                                               std::move(model),
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