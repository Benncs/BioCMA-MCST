#include <simulation/transport.hpp>
#include <simulation/update_flows.hpp>

namespace Simulation
{

  static void computeMatFlow(FlowInfo &flow, Simulation::MatFlow &matflow)
  {
    const auto mat_f_liq =
        Simulation::FlowmapToMat(flow.flows.data(), flow.flows.getN());
    const auto _mat_transition_liq =
        Simulation::get_transition_matrix(mat_f_liq);
    matflow.flows = mat_f_liq;
    matflow.transition_matrix = _mat_transition_liq;
  }

  static void computeMatFlow(std::span<double> flows,
                             size_t nc,
                             Simulation::MatFlow &matflow)
  {
    const auto mat_f_liq = Simulation::FlowmapToMat(flows, nc);
    const auto _mat_transition_liq =
        Simulation::get_transition_matrix(mat_f_liq);
    matflow.flows = mat_f_liq;
    matflow.transition_matrix = _mat_transition_liq;
  }

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   std::span<double> flows,
                   size_t nc,
                   Simulation::VecMatFlows &liq)
  {
    size_t current_index_mat = iteration_count % n_loop;
    auto &current_liq_matflow = liq.data[current_index_mat];

    if (iteration_count < n_loop)
    {
      computeMatFlow(flows, nc, current_liq_matflow);
    }

    unit.setLiquidFlow(&current_liq_matflow);

    iteration_count++;
  }

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   ReactorState *f,
                   Simulation::VecMatFlows &liq,
                   Simulation::VecMatFlows &gas)
  {
    size_t current_index_mat = iteration_count % n_loop;
    auto &current_liq_matflow = liq.data[current_index_mat];

    auto &current_gas_matflow = gas.data[current_index_mat];
    if (iteration_count < n_loop)
    {
      computeMatFlow(f->liquid_flow, current_liq_matflow);

      computeMatFlow(f->gas_flow, current_gas_matflow);
    }

    unit.mc_unit->domain.setLiquidNeighbors(f->liquid_flow.neigbors);
    unit.setLiquidFlow(&current_liq_matflow);

    unit.setGasFlow(&current_gas_matflow);

    iteration_count++;
  }

} // namespace Simulation