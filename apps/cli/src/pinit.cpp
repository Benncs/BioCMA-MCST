#include <pinit.hpp>

#include "flow_iterator.hpp"
#include "mpi.h"
#include <cstdio>
#include <simulation/models/simple_model.hpp>

#include <mc/mcinit.hpp>
#include <messages/pmessage.hpp>
#include <vector>

Simulation::SimulationUnit
pinit(ExecInfo &info, SimulationParameters &params, FlowIterator **it)
{

  std::vector<double> volumes;
  std::vector<std::vector<size_t>> liquid_neighbors;
  double d_t = 1e-1; // FIXME
  InitMessage payload;
  if (info.current_rank == 0)
  {
    auto ptr = new FlowIterator(params.flow_files);

    size_t n_t = params.final_time/d_t;

    ptr->setRepetition(n_t);

    *it = ptr;
    auto fstate = ptr->get();
    auto n_compartments = fstate->n_compartments;
    liquid_neighbors = ptr->get()->liquid_flow.neigbors;

    payload = {n_compartments, liquid_neighbors[0].size(), d_t};
    volumes = fstate->liquidVolume;
  }
  MPI_Bcast(&payload, sizeof(payload), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (info.current_rank != 0)
  {
    volumes.resize(payload.n_compartments);
    liquid_neighbors.resize(payload.n_compartments);
  }

  MPI_Bcast(volumes.data(), volumes.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // MPI_Bcast(neighbors.data(), neighbors.size(), MPI_INT, 0, MPI_COMM_WORLD);
  for (size_t i_n = 0; i_n < payload.n_compartments; ++i_n)
  {
    liquid_neighbors[i_n].resize(payload.n_neighbor);
    MPI_Bcast(liquid_neighbors[i_n].data(),
              liquid_neighbors[i_n].size(),
              MPI_INT,
              0,
              MPI_COMM_WORLD);
  }

  params.d_t = d_t;
  auto unit = MC::init_unit(info, volumes, std::move(liquid_neighbors));

  auto container = MC::init_container(info, params.n_particules);

  auto simulation = Simulation::SimulationUnit(params.n_species,
                                               std::move(unit),
                                               std::move(container),
                                               info,
                                               info.current_rank == 0);

  simulation.setVolumeLiquid(std::move(volumes));

  simulation.post_init(std::move(simple_model));

  return simulation;
}
