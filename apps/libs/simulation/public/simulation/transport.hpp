#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

#include <Eigen/Sparse>
#include <functional>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>

namespace Simulation
{
  Eigen::SparseMatrix<double>
  get_transition_matrix(const Eigen::SparseMatrix<double> &flows);


 
  struct MatFlow
  {
    Eigen::SparseMatrix<double> flows;
    Eigen::SparseMatrix<double> transition_matrix;

    MatFlow()=default;

    MatFlow(const Eigen::SparseMatrix<double> &_flows,
            const Eigen::SparseMatrix<double> &_tm)
    {
      flows = _flows;
      transition_matrix = _tm;
    }

    MatFlow(const Eigen::SparseMatrix<double> &&_flows,
            const Eigen::SparseMatrix<double> &&_tm)
    {
      flows = std::move(_flows);
      transition_matrix = std::move(_tm);
    }

    MatFlow &operator=(MatFlow &&rhs)
    {
      if (this != &rhs)
      {
        flows = std::move(rhs.flows);
        transition_matrix = std::move(rhs.transition_matrix);
      }
      return *this;
    }

    MatFlow(MatFlow &&other)
    {
      flows = std::move(other.flows);
      transition_matrix = std::move(other.transition_matrix);
    }
  };

  struct VecMatFlows
  {
    std::vector<MatFlow> data;

    VecMatFlows(size_t n)
    {
      data.resize(n);
    }
  };  


  using move_kernel =
      std::function<void(double, double, MC::Particles &, double)>;

  move_kernel pbf(MC::MonteCarloUnit &unit,
                  MC::ParticlesContainer &container,
                  const MatFlow &flows);
} // namespace Simulation

#endif //__TRANSPORT_HPP__