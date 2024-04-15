#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include <simulation/matflows.hpp>

namespace Simulation
{

  class ScalarSimulation
  {
  public:
    ScalarSimulation(ScalarSimulation &&other) noexcept;
    ScalarSimulation(size_t n_compartments, size_t n_species,size_t n_threads);
    Eigen::MatrixXd C;
    Eigen::MatrixXd Mtot;
    std::vector<Eigen::MatrixXd> contribs;

    Eigen::DiagonalMatrix<double, -1> V;

    void setV(std::vector<double> &&volumes);

    void merge(size_t i_thread);

    Eigen::MatrixXd biomass_contribution;
    inline size_t n_species() const
    {
      return this->n_r;
    }
    void performStep(double d_t,
                     const FlowMatrixType &m_transition,const Eigen::MatrixXd& transfer_g_l );

  private:
    void updateC();
    size_t n_r;
    size_t n_c;
  };

  inline void ScalarSimulation::setV(std::vector<double> &&volumes)
  {
    this->V.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));
  }

   inline void ScalarSimulation::merge(size_t i_thread)
  {
    this->biomass_contribution += this->contribs[i_thread];
  }

  inline ScalarSimulation* makeScalarSimulation(size_t n_compartments, size_t n_species,size_t n_threads)
  {
    return new ScalarSimulation(n_compartments, n_species,n_threads);
  }

} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__