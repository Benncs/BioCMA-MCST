#include <cma_read/flowmap.hpp>
#include <simulation/pc_hydro.hpp>
#include <simulation/update_flows.hpp>
#include <transport.hpp>

namespace Simulation
{
  static std::vector<double>
  compute_inverse_diagonal(std::span<const double> volumes)
  {
    std::vector<double> inverse_diagonal(volumes.size());

    std::transform(
        volumes.begin(),
        volumes.end(),
        inverse_diagonal.begin(),
        [](double volume)
        {
          if (volume == 0)
          {
            throw std::invalid_argument(
                "Setvolume: Null value of volume, matrix is not invertible");
          }
          return 1.0 / volume;
        });

    return inverse_diagonal;
  }

  // static void compute_MatFlow(const FlowInfo &flow,
  // Simulation::PreCalculatedHydroState &matflow)
  // {
  //   auto view = flow.getViewFlows();

  //   const auto _mat_transition_liq =
  //       Simulation::get_transition_matrix(view);

  //   matflow.transition_matrix = _mat_transition_liq;
  // }

  static void compute_MatFlow(const FlowMap::FlowMap_const_view_t &flows_view,
                              Simulation::PreCalculatedHydroState &matflow)
  {
    const auto _mat_transition_liq =
        Simulation::get_transition_matrix(flows_view);
    matflow.transition_matrix = _mat_transition_liq;
  }

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   std::span<double> flows,
                   size_t nc,
                   Simulation::BasicCacheHydro &liq)
  {
    size_t current_index_mat = iteration_count % n_loop;
    auto &current_liq_matflow = liq.data[current_index_mat];
    const auto mat_f_liq_view = FlowMap::FlowMap_const_view_t(flows, nc);
    if (iteration_count < n_loop)
    {
      compute_MatFlow(mat_f_liq_view, current_liq_matflow);
      current_liq_matflow.cumulative_probability =
          get_CP(unit.mc_unit->domain.getNeighbors(),
                 current_liq_matflow.transition_matrix);
    }

    unit.setLiquidFlow(&current_liq_matflow);

    iteration_count++;
  }

  void update_flow(size_t &iteration_count,
                   size_t n_per_flowmap,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   const ReactorState &reactor_state,
                   Simulation::BasicCacheHydro &liquid_flows,
                   Simulation::BasicCacheHydro &gas_flows,
                   bool tpf)
  {
    static size_t element_count = 0;

    size_t current_index_mat = iteration_count % n_loop;
    auto &current_liq_matflow = liquid_flows.data[current_index_mat];
    unit.mc_unit->domain.setLiquidNeighbors(
        reactor_state.liquid_flow.getViewNeighors());
    auto &current_gas_matflow = gas_flows.data[current_index_mat];
    if (iteration_count < n_loop)
    {
      compute_MatFlow(reactor_state.liquid_flow.getViewFlows(),
                      current_liq_matflow);

      current_liq_matflow.inverse_volume =
          compute_inverse_diagonal(reactor_state.liquidVolume);
      if (tpf)
      {
        compute_MatFlow(reactor_state.gas_flow.getViewFlows(),
                        current_gas_matflow);
        current_gas_matflow.inverse_volume =
            compute_inverse_diagonal(reactor_state.gasVolume);
      }

      current_liq_matflow.cumulative_probability =
          get_CP(unit.mc_unit->domain.getNeighbors(),
                 current_liq_matflow.transition_matrix);
    }

    unit.setLiquidFlow(&current_liq_matflow);
    unit.setGasFlow(&current_gas_matflow);

    if (++element_count == n_per_flowmap)
    {
      iteration_count++;
      element_count = 0;
    }
  }

} // namespace Simulation