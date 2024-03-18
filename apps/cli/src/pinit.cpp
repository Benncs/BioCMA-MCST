#include <pinit.hpp>

#include "flow_iterator.hpp"
#include "mpi.h"
#include <cstdio>

#include <mc/mcinit.hpp>
#include <messages/pmessage.hpp>
#include <vector>

Simulation::SimulationUnit
pinit(ExecInfo &info, SimulationParameters &params, FlowIterator **it)
{

  std::vector<double> volumes;
  std::vector<std::vector<size_t>> neighbors;
  double d_t = 1e-1; // FIXME
  InitMessage payload;
  if (info.current_rank == 0)
  {
    auto ptr = new FlowIterator(params.flow_files);

    ptr->setRepetition(500);

    *it = ptr;
    auto n_compartments = ptr->get()->flows.getN();
    neighbors = ptr->get()->neigbors;

    payload = {n_compartments, neighbors.size() / n_compartments, d_t};
    volumes = std::vector(n_compartments, 0.2);
  }
  MPI_Bcast(&payload, sizeof(payload), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (info.current_rank != 0)
  {
    volumes.resize(payload.n_compartments);
    neighbors.resize(payload.n_compartments);
  }

  MPI_Bcast(volumes.data(), volumes.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // MPI_Bcast(neighbors.data(), neighbors.size(), MPI_INT, 0, MPI_COMM_WORLD);
  for (size_t i_n = 0; i_n < payload.n_compartments; ++i_n)
  {
    neighbors[i_n].resize(payload.n_neighbor);
    MPI_Bcast(neighbors[i_n].data(),
              neighbors[i_n].size(),
              MPI_INT,
              0,
              MPI_COMM_WORLD);
  }

  params.d_t = d_t;
  auto unit = MC::init_unit(info, volumes, std::move(neighbors));

  auto container = MC::init_container(info, params.n_particules);

  auto simulation = Simulation::SimulationUnit(params.n_species,
                                               std::move(unit),
                                               std::move(container),
                                               info.current_rank == 0);

  return simulation;
}
