#include <simulation/transport.hpp>
#include <simulation/update_flows.hpp>

namespace Simulation
{
  static std::vector<double> compute_inverse_diagonal(std::span<double> volumes) {
    std::vector<double> inverse_diagonal(volumes.size());
    
    std::transform(volumes.begin(), volumes.end(), inverse_diagonal.begin(), [](double volume) {
        if (volume == 0) {
            throw std::invalid_argument("Setvolume: Null value of volume, matrix is not invertible");
        }
        return 1.0 / volume;
    });

    return inverse_diagonal;
}

  static void compute_MatFlow(FlowInfo &flow, Simulation::MatFlow &matflow)
  {
    const auto mat_f_liq =
        Simulation::FlowmapToMat(flow.flows.data(), flow.flows.getN());
    const auto _mat_transition_liq =
        Simulation::get_transition_matrix(mat_f_liq);




    matflow.flows = mat_f_liq;
    matflow.transition_matrix = _mat_transition_liq;
  }

  static void compute_MatFlow(std::span<double> flows,
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
                   Simulation::BasicCacheMatflows &liq)
  {
    size_t current_index_mat = iteration_count % n_loop;
    auto &current_liq_matflow = liq.data[current_index_mat];

    if (iteration_count < n_loop)
    {
      compute_MatFlow(flows, nc, current_liq_matflow);
    }

    unit.setLiquidFlow(&current_liq_matflow);

    iteration_count++;
  }

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   ReactorState &reactor_state,
                   Simulation::BasicCacheMatflows &liquid_flows,
                   Simulation::BasicCacheMatflows &gas_flows)
  {
    size_t current_index_mat = iteration_count % n_loop;
    auto &current_liq_matflow = liquid_flows.data[current_index_mat];

    auto &current_gas_matflow = gas_flows.data[current_index_mat];
    if (iteration_count < n_loop)
    {
      compute_MatFlow(reactor_state.liquid_flow, current_liq_matflow);
     
      compute_MatFlow(reactor_state.gas_flow, current_gas_matflow);
  
      current_liq_matflow.inverse_volume = compute_inverse_diagonal(reactor_state.liquidVolume);
      current_gas_matflow.inverse_volume = compute_inverse_diagonal(reactor_state.gasVolume);
    }

    unit.mc_unit->domain.setLiquidNeighbors(std::cref(reactor_state.liquid_flow.neigbors));
    unit.setLiquidFlow(&current_liq_matflow);
    unit.setGasFlow(&current_gas_matflow);
   

    iteration_count++;
  }

} // namespace Simulation