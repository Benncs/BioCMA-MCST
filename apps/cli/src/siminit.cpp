

#include <cstddef>
#include <cstdio>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cma_read/flow_iterator.hpp>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/neighbors.hpp>
#include <cmt_common/cma_case.hpp>
#include <common/common.hpp>
#include <mc/mcinit.hpp>
#include <mc/prng/distribution.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <simulation/simulation.hpp>

#include <simulation/update_flows.hpp>

static CmaRead::ReactorState const *
init_state(SimulationParameters &params,
           std::unique_ptr<CmaRead::FlowIterator> &flow_handle,
           const CmtCommons::cma_exported_paths_t &paths);

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
                std::unique_ptr<Simulation::FlowMapTransitioner> &transitioner,
                KModel model)
{

  const auto &user_params = params.user_params;

  std::vector<double> liq_volume;
  std::vector<double> gas_volume;
  std::vector<size_t> worker_neighbor_data;
  CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors;
  std::unique_ptr<CmaRead::FlowIterator> _flow_handle = nullptr;
  init_host_only(info,
                 params,
                 _flow_handle,
                 liquid_neighbors,
                 liq_volume,
                 gas_volume,
                 worker_neighbor_data);
  if constexpr (RT::use_mpi)
  {
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
  }

  // MPI_W::bcst_iterator(_flow_handle, info.current_rank);
  if (info.current_rank != 0)
  {

    auto n_col = worker_neighbor_data.size() / params.n_compartments;
    liquid_neighbors = CmaRead::Neighbors::Neighbors_const_view_t(
        worker_neighbor_data, params.n_compartments, n_col);
  }
  liquid_neighbors.set_row_major();

  MC::UniformLawINT law_param = {0, static_cast<int>(params.n_compartments-1)}; 

  // MPI_W::barrier(); //Useless ?
  auto mc_unit = MC::init(model.init_kernel,
                          info,
                          user_params.numper_particle,
                          liq_volume,
                          liquid_neighbors,
                          law_param);

  bool f_init_gas_flow = info.current_rank == 0 && params.is_two_phase_flow;

  auto simulation =
      std::make_unique<Simulation::SimulationUnit>(info,
                                                   std::move(mc_unit),
                                                   gas_volume,
                                                   liq_volume,
                                                   params.n_species,
                                                   model,
                                                   f_init_gas_flow);

  // Calculate the total number of time steps
  const auto n_t = static_cast<size_t>(user_params.final_time / params.d_t) + 1;

  transitioner = std::make_unique<Simulation::FlowMapTransitioner>(
      params.n_different_maps,
      params.n_per_flowmap,
      Simulation::FlowMapTransitioner::Discontinuous,
      n_t,
      std::move(_flow_handle),
      f_init_gas_flow);
  // params.is_two_phase_flow);

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
  std::string case_name = params.user_params.cma_case_path + "/cma_case";
  CmtCommons::CMACaseInfo cma_case =
      CmtCommons::CMACaseInfoReader::load_case(case_name);

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
    params.d_t = (params.user_params.delta_time == 0.)
                     ? _flow_handle->MinLiquidResidenceTime() / 100.
                     : params.user_params.delta_time;
  }
  else
  {
    params.d_t = params.user_params.delta_time;
  }

  // const auto n_t = static_cast<size_t>(params.final_time / params.d_t) + 1;

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

  // _flow_handle->setRepetition(n_repetition, n_per_flowmap);
  register_run(info, params);
}

static CmaRead::ReactorState const *
init_state(SimulationParameters &params,
           std::unique_ptr<CmaRead::FlowIterator> &flow_handle,
           const CmtCommons::cma_exported_paths_t &paths)
{
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
  }
  catch (const std::exception &e)
  {
    std::stringstream err;
    err << "Error while reading files\t:" << e.what();

    throw std::runtime_error(err.str());
  }

  if (state == nullptr)
  {
    throw std::runtime_error("Error while reading files");
  }
  return state;
}
