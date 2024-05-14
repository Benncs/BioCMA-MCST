#ifndef __SCALAR_SIMULATION_HPP__
#define __SCALAR_SIMULATION_HPP__

#include <simulation/pc_hydro.hpp>
#include <span>

//TODO REMOVE 
#include <iostream>

namespace Simulation
{

  class ScalarSimulation
  {
  public:
    ScalarSimulation(ScalarSimulation &&other) noexcept;
    ScalarSimulation(const ScalarSimulation &other) noexcept = delete;
    ScalarSimulation(size_t n_compartments, size_t n_species, size_t n_threads,std::span<double> volume);
    ScalarSimulation operator=(const ScalarSimulation &other) = delete;
    ScalarSimulation operator=(ScalarSimulation &&other) = delete;

    ~ScalarSimulation() = default;

    Eigen::MatrixXd concentration;
    Eigen::MatrixXd total_mass;

    // Getters
    [[nodiscard]] std::span<double const> getVolumeData() const;

    [[nodiscard]] std::span<double> getContributionData();

    std::span<Eigen::MatrixXd> getThreadContribs();

    Eigen::DiagonalMatrix<double, -1> &getVolume();

    std::span<double> getCData();

    [[nodiscard]] inline size_t n_species() const;

    // Setters

    void setVolumes(std::span<double> volumes, std::span<double> inv_volumes);

    void merge(size_t i_thread);

    Eigen::ArrayXXd vec_kla; // TODO : Clean this

    Eigen::MatrixXd biomass_contribution;

    void performStep(double d_t,
                     const FlowMatrixType &m_transition,
                     const Eigen::MatrixXd &transfer_gas_liquid);

  private:
    Eigen::DiagonalMatrix<double, -1> volumes_inverse;
    std::vector<Eigen::MatrixXd> contribs;
    Eigen::DiagonalMatrix<double, -1> m_volumes;

    void updateC();
    size_t n_r;
    size_t n_c;
  };

  inline std::span<Eigen::MatrixXd> ScalarSimulation::getThreadContribs()
  {
    return contribs;
  }

  inline Eigen::DiagonalMatrix<double, -1> &ScalarSimulation::getVolume()
  {
    return m_volumes;
  }

  inline size_t ScalarSimulation::n_species() const
  {
    return this->n_r;
  }

  inline std::span<double> ScalarSimulation::getCData()
  {

    return {this->concentration.data(), static_cast<size_t>(this->concentration.size())};
  }

  inline std::span<double> ScalarSimulation::getContributionData()
  {
    return {this->biomass_contribution.data(),
            static_cast<size_t>(this->biomass_contribution.size())};
  }

  inline std::span<double const> ScalarSimulation::getVolumeData() const
  {
    return {m_volumes.diagonal().data(), static_cast<size_t>(m_volumes.rows())};
  }


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
                                                size_t n_threads,std::span<double> volumes)
  {
    return new ScalarSimulation(n_compartments, n_species, n_threads,volumes);
  }

} // namespace Simulation

#endif //__SCALAR_SIMULATION_HPP__