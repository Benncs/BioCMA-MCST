#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include "cma_read/light_2d_view.hpp"
#include "mc/domain.hpp"
#include <Kokkos_Core.hpp>
#include <cma_read/flowmap.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <simulation/pc_hydro.hpp>

namespace Simulation
{

  FlowMatrixType
  get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows);

  std::vector<double> get_diag_transition(const FlowMatrixType &m_transition);

  Eigen::MatrixXd get_CP(CmaRead::Neighbors::Neighbors_const_view_t neighbors,
                         const FlowMatrixType &m_transition);

  void kernel_exit(double d_t,
                   double random_number,
                   Kokkos::View<MC::ContainerState *, Kokkos::LayoutStride> domain,
                   MC::Particles &particle);

  void kernel_move(
      double random_number,
      double random_number2,
      MC::Particles &particle,
      double d_t,
      std::span<const double> diag_transition,
      CmaRead::L2DView<const double> cumulative_probability,
      Kokkos::View<MC::ContainerState *, Kokkos::LayoutStride> domain_view,
      Kokkos::View<const size_t **, Kokkos::LayoutStride> view_neighbors);

} // namespace Simulation

#endif //__TRANSPORT_HPP__