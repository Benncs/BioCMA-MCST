#include "messages/impl_op.hpp"
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <host_specific.hpp>
#include <messages/wrap_mpi.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <sync.hpp>

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
               std::shared_ptr<FlowIterator> _flow_handle)
{

  size_t iteration_count = 0;
  size_t n_loop = params.n_different_maps;

  auto liquid_flows = Simulation::VecMatFlows(n_loop);
  auto gas_flows = Simulation::VecMatFlows(n_loop);

  double d_t = params.d_t;
    std::cout<<params.final_time<<" "<<d_t<<std::endl;
  for (auto &&f : *_flow_handle)
  {
#ifdef DEBUG
    static size_t __debug_loop_cntr = 0;
    std::cout << __debug_loop_cntr << std::endl;
    __debug_loop_cntr++;
#endif

    simulation.state = &f;

    host_dispatch(exec, MPI_W::SIGNALS::RUN, f.liquid_flow.flows.data());

    // Send the size of the neighbor vectors to all processes
    size_t neighbor_size = f.liquid_flow.neigbors[0].size();

    // Send each neighbor vector to all processes
    for (const auto &neighbor : f.liquid_flow.neigbors)
    {
      const unsigned long *buf = neighbor.data();
      int size = neighbor.size();
      for (int j = 1; j < static_cast<int>(exec.n_rank); ++j)
      {
        // MPI_W::send(std::move(neighbor_size), j);
        // MPI_W::send_v<size_t>(neighbor, j);
        MPI_Send(&neighbor_size, 1, MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD);
        MPI_Send(buf, size, MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD);
      }
    }

    simulation.mc_unit->domain.setLiquidNeighbors(f.liquid_flow.neigbors);

    Simulation::update_flow(
        iteration_count, n_loop, simulation, &f, liquid_flows, gas_flows);

    simulation.cycleProcess(d_t);

    sync_step(exec, simulation);
    simulation.step(d_t);
    sync_prepare_next(exec, simulation);

  }
}
