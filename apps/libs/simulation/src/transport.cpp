#include "common/thread_safe_container.hpp"
#include "mc/domain.hpp"
#include <simulation/pc_hydro.hpp>
#include <get_cumulative_proba.hpp>
#include <transport.hpp>

// TODO REMOVE
#include <iostream>

// TODO REMOVE
template <typename T>
using RowMajorDynMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using ColMajorDynMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

namespace Simulation
{

  static size_t
  find_next_compartment(int i_compartment,
                        double random_number,
                        std::span<const size_t> i_neighbor,
                        const Eigen::MatrixXd &cumulative_probability);

  // move_kernel_t population_balance_flow(MC::ReactorDomain &domain,

  //                                       const PreCalculatedHydroState *flows)
  // {

  //   const auto cumulative_probability = get_CP(domain.getNeighbors(), *flows);

  //   const auto &m_transition = flows->transition_matrix;

  //   auto move_kernel = [m_transition,
  //                       cumulative_probability](double random_number,
  //                                               double random_number2,
  //                                               MC::ReactorDomain &domain,
  //                                               MC::Particles &particle,
  //                                               double d_t) -> void
  //   {
  //     kernel_move(random_number,
  //                 random_number2,
  //                 domain,
  //                 particle,
  //                 d_t,
  //                 m_transition,
  //                 cumulative_probability);
  //   };

  //   return move_kernel;
  // }

  Eigen::MatrixXd get_CP(CmaRead::Neighbors::Neighbors_const_view_t neighbors,
                         const FlowMatrixType &m_transition)
  {
    size_t n_zone = neighbors.getNRow();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(
        static_cast<int>(n_zone), static_cast<int>(neighbors.getNCol()));

    for (int k_compartment = 0; k_compartment < static_cast<int>(n_zone);
         ++k_compartment)
    {
      double cumsum = 0.;
      int count_nei = 0;
      auto row = neighbors.getRow(k_compartment);
      for (auto &&i_neighbor : row)
      {
        // std::cout<<i_neighbor<<"\t";
        auto colId = static_cast<int>(i_neighbor);
        if (colId == k_compartment)
        {
          break;
        }

        const double proba_out =
            m_transition.coeff(k_compartment, colId) /
            std::abs(m_transition.coeff(k_compartment, k_compartment));

        const double p_cp = proba_out + cumsum;
        P.coeffRef(k_compartment, count_nei) = p_cp;
        count_nei++;
        cumsum += proba_out;
      }
      // std::cout<<"\r\n";
    }
    return P;
  }

  void kernel_move(double random_number,
                   double random_number2,
                   MC::ReactorDomain &domain,
                   MC::Particles &particle,
                   double d_t,
                   const FlowMatrixType &m_transition,
                   const Eigen::MatrixXd &cumulative_probability)
  {
    const size_t i_compartment = particle.current_container;
    const int rowId = static_cast<int>(i_compartment);

    auto &current_container = domain[i_compartment];
    const std::span<const size_t> i_neighbor =
        domain.getNeighbors(i_compartment);
    
    const double &v_p = current_container.volume_liq;

    const double leaving_flow = std::abs(m_transition.coeff(rowId, rowId));

    const double theta_p = v_p / leaving_flow;

    const double probability = 1 - std::exp(-d_t / theta_p);

    if (random_number >= probability)
    {
      return;
    }

    size_t next = find_next_compartment(
        rowId, random_number2, i_neighbor, cumulative_probability);

    __ATOM_DECR__(current_container.n_cells) // Cell leaves current
                                             // compartment
    __ATOM_INCR__(domain[next].n_cells);     // Cell go to new compartment
    particle.current_container = next;
  }

  size_t find_next_compartment(int i_compartment,
                               double random_number,
                               std::span<const size_t> i_neighbor,
                               const Eigen::MatrixXd &cumulative_probability)
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor[0];
   
    for (int k_neighbor = 0; k_neighbor < max_neighbor - 1; ++k_neighbor)
    {
      auto pi = cumulative_probability.coeff(i_compartment, k_neighbor);
      auto pn = cumulative_probability.coeff(i_compartment, k_neighbor + 1);
      if (random_number <= pn && pi <= random_number)
      {
        next = i_neighbor[k_neighbor + 1];

        // break;
      }
    }

    return next;
  }

  
   FlowMatrixType get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows)
  {
    int n_compartiments = static_cast<int>(flows.getNRow()); // It SHOULD be square
    auto rd = flows.data();
        // auto rd = flows.data();
    std::vector<double> flow_copy = std::vector<double>(rd.begin(),rd.end());
    FlowMatrixType m_transition = Eigen::Map<Eigen::MatrixXd>(flow_copy.data(),n_compartiments,n_compartiments);
    // FlowMatrixType m_transition(n_compartiments, n_compartiments);
    

    for(int i =0;i<n_compartiments;++i)
    {
      m_transition.coeffRef(i,i) = -1.*m_transition.col(i).sum();
    }

    // TODO OPTI
    // for (int i = 0; i < n_compartiments; ++i)
    // {
    //   double sum = 0.;
    //   for (int j = 0; j < n_compartiments; ++j)
    //   {
    //     sum += flows(i, j);
    //     if (i != j)
    //     {
    //       m_transition.coeffRef(i, j) = flows(i, j);
    //     }
    //   }
    //   m_transition.coeffRef(i, i) = -sum;
    // }

    return m_transition;
  }

} // namespace Simulation
