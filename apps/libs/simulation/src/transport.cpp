#include "common/thread_safe_container.hpp"
#include "mc/domain.hpp"
#include "mc/particles/mcparticles.hpp"
#include <cmt_common/zip.hpp>
#include <get_cumulative_proba.hpp>
#include <simulation/pc_hydro.hpp>
#include <stdexcept>
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
                        CmaRead::L2DView<const double> cumulative_probability);

  // move_kernel_t population_balance_flow(MC::ReactorDomain &domain,

  //                                       const PreCalculatedHydroState *flows)
  // {

  //   const auto cumulative_probability = get_CP(domain.getNeighbors(),
  //   *flows);

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
    }
    return P;
  }
  inline bool probability_leaving(double random_number,
                                  double volume,
                                  double flow,
                                  double dt)
  {
    // const double theta_p = volume / flow;
    // const double probability = 1 - std::exp(-dt / theta_p);
    // bool p1 = random_number < probability;

    // return (dt/theta_p)>(-std::log(1-random_number));
    return (dt*flow/volume)>(-std::log(1-random_number));

    // if(p1!=p2)
    // {
    //   throw std::runtime_error("ALED");
    // }

    // return p2;
  }

  static const std::vector<size_t> v_tmp_index_leaving_flow = {10};
  static const std::vector<double> v_tmp_leaving_flow = {0.000011758 / 10};

  static auto [zip_it,zip_end] = CmtCommons::zip(v_tmp_index_leaving_flow,v_tmp_leaving_flow);

  void kernel_exit(double d_t,
                   double random_number,
                   MC::ReactorDomain &domain,
                   MC::Particles &particle)
  {
    
    
    const auto lambda = [&particle, random_number, &domain, d_t](auto &&index, auto &&flow)
        {
          if (particle.current_container != index ||
              particle.status != MC::CellStatus::IDLE)
          {
            return;
          }

          if (probability_leaving(
                  random_number, domain[index].volume_liq, flow, d_t))
          {
            particle.status = MC::CellStatus::OUT;
          }
        };

    for(size_t i =0;i<v_tmp_index_leaving_flow.size();++i)
    {
      auto& index = v_tmp_index_leaving_flow[i];
      auto& flow = v_tmp_leaving_flow[i];
      lambda(index,flow);
    }

    // CmtCommons::foreach_zip(
    //     lambda,
    //     v_tmp_index_leaving_flow,
    //     v_tmp_leaving_flow);
  }

  void kernel_move(double random_number,
                   double random_number2,
                   MC::ReactorDomain &domain,
                   MC::Particles &particle,
                   double d_t,
                   std::span<const double> diag_transition,
                   CmaRead::L2DView<const double> cumulative_probability)
  {
    const size_t i_compartment = particle.current_container;
    const int rowId = static_cast<int>(i_compartment);

    auto &current_container = domain[i_compartment];
    const std::span<const size_t> i_neighbor =
        domain.getNeighbors(i_compartment);

    if (!probability_leaving(random_number,
                             current_container.volume_liq,
                             diag_transition[i_compartment],
                             d_t))
    {
      return;
    }

    const size_t next = find_next_compartment(
        rowId, random_number2, i_neighbor, cumulative_probability);

    __ATOM_DECR__(current_container.n_cells) // Cell leaves current
                                             // compartment
    __ATOM_INCR__(domain[next].n_cells);     // Cell go to new compartment
    particle.current_container = next;
  }

  size_t
  find_next_compartment(int i_compartment,
                        double random_number,
                        std::span<const size_t> i_neighbor,
                        CmaRead::L2DView<const double> cumulative_probability)
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor[0];

    for (int k_neighbor = 0; k_neighbor < max_neighbor - 1; ++k_neighbor)
    {
      const auto pi = cumulative_probability(i_compartment, k_neighbor);
      const auto pn = cumulative_probability(i_compartment, k_neighbor + 1);
      if (random_number <= pn && pi <= random_number)
      {
        next = i_neighbor[k_neighbor + 1];

        // break;
      }
    }

    return next;
  }

  std::vector<double> get_diag_transition(const FlowMatrixType &m_transition)
  {
    std::vector<double> res(m_transition.rows());
    for (int i = 0; i < m_transition.rows(); ++i)
    {
      res[i] = -m_transition.coeff(i, i);
    }
    return res;
  }

  FlowMatrixType
  get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows)
  {
    int n_compartments =
        static_cast<int>(flows.getNRow()); // It SHOULD be square

    // Uncomment with dense matrix
    //  auto rd = flows.data();
    //  std::vector<double> flow_copy =
    //  std::vector<double>(rd.begin(),rd.end()); FlowMatrixType m_transition =
    //  Eigen::Map<Eigen::MatrixXd>(flow_copy.data(),n_compartments,n_compartments);

    FlowMatrixType m_transition =
        FlowMatrixType(n_compartments, n_compartments);

    for (size_t i = 0; i < n_compartments; ++i)
    {
      for (size_t j = 0; j < n_compartments; ++j)
      {
        const auto val = flows(i, j);
        if (val != 0.)
        {
          m_transition.coeffRef(i, j) += val;
        }
      }
    }

    for (int i = 0; i < n_compartments; ++i)
    {
      m_transition.coeffRef(i, i) = -1. * m_transition.col(i).sum();
    }

    m_transition.makeCompressed();

    return m_transition;
  }

} // namespace Simulation
