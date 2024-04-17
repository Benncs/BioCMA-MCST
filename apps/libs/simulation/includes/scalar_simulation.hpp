#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include <simulation/matflows.hpp>
#include <span>

namespace Simulation
{

  class ScalarSimulation
  {
  public:
    ScalarSimulation(ScalarSimulation &&other) noexcept;
    ScalarSimulation(const ScalarSimulation &other) noexcept = delete;
    ScalarSimulation(size_t n_compartments, size_t n_species, size_t n_threads);
    ScalarSimulation operator=(const ScalarSimulation& other) =delete;
    ScalarSimulation operator=(ScalarSimulation&& other) =delete;

    ~ScalarSimulation() = default;

    Eigen::MatrixXd C;
    Eigen::MatrixXd Mtot;
    std::vector<Eigen::MatrixXd> contribs;

    Eigen::DiagonalMatrix<double, -1> m_volumes;

    void setVolumes(std::span<double> volumes, std::span<double> inv_volumes);

    void merge(size_t i_thread);

    Eigen::ArrayXXd vec_kla; // TODO : Clean this

    Eigen::MatrixXd biomass_contribution;
    [[nodiscard]] inline size_t n_species() const
    {
      return this->n_r;
    }
    void performStep(double d_t,
                     const FlowMatrixType &m_transition,
                     const Eigen::MatrixXd &transfer_gas_liquid);

  private:
    Eigen::DiagonalMatrix<double, -1> volumes_inverse;

    void updateC();
    size_t n_r;
    size_t n_c;
  };

  inline void ScalarSimulation::setVolumes(std::span<double> volumes,
                                           std::span<double> inv_volumes)
  {

    // SIGFAULT ?
    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    this->volumes_inverse.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        inv_volumes.data(), static_cast<int>(inv_volumes.size()));
  }

  inline void ScalarSimulation::merge(size_t i_thread)
  {
    this->biomass_contribution += this->contribs[i_thread];
  }

  inline ScalarSimulation *makeScalarSimulation(size_t n_compartments,
                                                size_t n_species,
                                                size_t n_threads)
  {
    return new ScalarSimulation(n_compartments, n_species, n_threads);
  }

} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__