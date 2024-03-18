#include <get_cumulative_proba.hpp>
#include <random>
#include <simulation/transport.hpp>
  #include <iostream>

namespace Simulation
{

  Eigen::SparseMatrix<double>
  get_transition_matrix(const Eigen::SparseMatrix<double> &flows)
  {
    auto n_compartiments = flows.rows(); // It's SHOULD be square
    Eigen::SparseMatrix<double> m_transition(n_compartiments, n_compartiments);

    // TODO OPTI
    for (int i = 0; i < n_compartiments; ++i)
    {
      for (int j = 0; j < n_compartiments; ++j)
      {
        if (i != j)
        {
          m_transition.coeffRef(i, j) = flows.coeff(i, j);
        }
        else
        {
          Eigen::RowVectorXd row = flows.row(i);

          auto sum_i =
              std::accumulate(row.begin(), row.end(), 0.) - flows.coeff(i, i);

          m_transition.coeffRef(i, i) = -sum_i;
        }
      }
    }
    m_transition.makeCompressed();
    return m_transition;
  }

  // def pbf(domain,container,f,m,dt:float):

  //     N_NEIGHBOR = 6  # CUBE
  //     n_c = domain.size
  //     pik = f/np.sum(f)
  //     CP = np.zeros((n_c,N_NEIGHBOR))
  //     cs = np.cumsum(pik,axis=1)
  //     for i in range(N_NEIGHBOR):
  //         CP[:,i] = cs[:,i]

  //     f_move(container.to_process._data,domain,m,dt,CP)
  move_kernel pbf(MC::MonteCarloUnit &unit,
                  MC::ParticlesContainer &container,
                  const MatFlow &flows)
  {
    auto n_c = static_cast<int>(unit.domain.n_compartments());


    auto cumulative_propability =
        get_CP(unit.domain.getNeighbors(), n_c, flows);

    auto &m_transition = flows.transition_matrix;

    auto move_kernel =
        [&unit, &m_transition, &cumulative_propability](double random_number,
                                                        double random_number2,
                                                        MC::Particles &particle,
                                                        double d_t) -> void
    {
      const size_t i_compartment = particle.current_container;
      auto &current_container = unit.domain[i_compartment];
      const auto i_neighbor = unit.domain.getNeighbors();
      const int max_neighbor = static_cast<int>(i_neighbor[0].size());

      const double &v_p = current_container.volume;
      auto c = current_container.concentrations;

      
      const double leaving_flow = m_transition.coeff(
          static_cast<int>(i_compartment), static_cast<int>(i_compartment));

      double theta_p = -v_p / leaving_flow;

      double probability = std::exp(-d_t / theta_p);

      bool leaving = false;
      
      if (random_number > probability)
      {
        return;
      }
      else
      {
        leaving = true;
        size_t next = i_neighbor[i_compartment][0];
        for (int k = 1; k < max_neighbor - 1; ++k)
        {
          if (cumulative_propability.coeff(static_cast<int>(i_compartment), k) <
                  random_number2 &&
              random_number2 < cumulative_propability.coeff(
                                   static_cast<int>(i_compartment), k + 1))
          {
            next = i_neighbor[i_compartment][k];
          }
        }
#pragma omp critical
        {
          current_container.n_cells -= static_cast<size_t>(leaving);
          particle.current_container = next;
          unit.domain[next].n_cells += 1;
        }
      }
    };

    return move_kernel;
  }

} // namespace Simulation
