#include "common/thread_safe_container.hpp"
#include "mc/prng/prng.hpp"
#include "models/types.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <mc/domain.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/thread_private_data.hpp>
#include <mc/unit.hpp>

#include <memory>

#include <scalar_simulation.hpp>
#include <simulation/simulation.hpp>
#include <stdexcept>
#include <transport.hpp>

#ifndef USE_PYTHON_MODULE
#  include <omp.h>
#else
#  define omp_get_thread_num() 0 // 1 thread
#endif

// TODO REMOVE
#include <iostream>
#include <utility>

namespace Simulation
{
  void pimpl_deleter::operator()(ScalarSimulation *ptr) const
  {
    delete ptr;
  }

  void initF(pimp_ptr_t &liq, pimp_ptr_t &gas)
  {
    auto rng = MC::PRNG();
    bool host = gas != nullptr;
    // for (auto i = 0; i < liq->concentration.cols(); ++i)
    // {
    //   // liq->concentration.coeffRef(1, i) = 0;
    //   // liq->concentration.coeffRef(0, static_cast<int>(i)) = 0.5; // 0.5
    //   G/LP
    //   // Glucose
    //   double rngn = (i < liq->concentration.cols() / 3)
    //                     ? rng.uniform_double_rand(0, 5)
    //                     : rng.uniform_double_rand(0, 20);
    //   liq->concentration.coeffRef(0, static_cast<int>(i)) = rngn;
    //   if (host)
    //   {
    //     gas->concentration.coeffRef(1, i) =
    //         15e-3 * 1e5 / (8.314 * (273.15 + 30)) * 0.2;
    //   }
    // }

    liq->concentration.coeffRef(0, static_cast<int>(0)) =
        100000; // 0.5 G/L Glucose

    liq->total_mass = liq->concentration * liq->getVolume();
    if (host)
    {
      gas->total_mass = gas->concentration * gas->getVolume();
    }
  }

  void p_kernel(double d_t,
                MC::MonteCarloUnit &unit,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::ThreadPrivateData> _extras,
                const KModel &_kmodel,
                MC::Particles &p,
                auto &m_transition,
                auto &cumulative_probability);

  SimulationUnit::SimulationUnit(SimulationUnit &&other) noexcept
      : mc_unit(std::move(other.mc_unit)),
        is_two_phase_flow(other.is_two_phase_flow), n_thread(other.n_thread),
        flow_liquid(other.flow_liquid), flow_gas(other.flow_gas),
        kmodel(other.kmodel)
  {
  }

  SimulationUnit::SimulationUnit(
      const ExecInfo &info,
      std::unique_ptr<MC::MonteCarloUnit> &&_unit,
      std::span<double> volumesgas,
      std::span<double> volumesliq,
      size_t n_species,
      KModel _km,
      MC::DistributionVariantInt &&distribution_variant,
      bool _gas_flow)
      : mc_unit(std::move(_unit)), is_two_phase_flow(_gas_flow),
        n_thread(info.thread_per_process), flow_liquid(nullptr),
        flow_gas(nullptr), kmodel(_km)
  {
    if (this->mc_unit->extras.empty())
    {
      throw std::runtime_error("Extra particle not initialised");
    }

    this->liquid_scalar = std::unique_ptr<ScalarSimulation, pimpl_deleter>(
        makeScalarSimulation(mc_unit->domain.getNumberCompartments(),
                             n_species,
                             n_thread,
                             volumesliq));

    this->gas_scalar =
        (is_two_phase_flow)
            ? std::unique_ptr<ScalarSimulation, pimpl_deleter>(
                  makeScalarSimulation(mc_unit->domain.getNumberCompartments(),
                                       n_species,
                                       0,
                                       volumesgas)) // No contribs for gas
            : nullptr;

    // FIXME
    initF(liquid_scalar, gas_scalar);

    post_init_container(std::move(distribution_variant));
    post_init_compartments();
  }

  void SimulationUnit::setVolumes(std::span<double> volumesgas,
                                  std::span<double> volumesliq)
  {

    std::span<double const> vg;
    this->liquid_scalar->setVolumes(volumesliq, flow_liquid->inverse_volume);
    if (gas_scalar)
    {
      this->gas_scalar->setVolumes(volumesgas, flow_gas->inverse_volume);
      vg = gas_scalar->getVolumeData();
    }
    else
    {
      vg = volumesgas;
    }

    std::span<double const> vl = liquid_scalar->getVolumeData();

    this->mc_unit->domain.setVolumes(vg, vl);
  }

  void SimulationUnit::post_init_compartments()
  {
    int i = 0;
    size_t n_species = this->liquid_scalar->n_species();
    std::for_each(mc_unit->domain.begin(),
                  mc_unit->domain.end(),
                  [this, n_species, &i](auto &&cs)
                  {
                    cs.concentrations = std::span(
                        this->liquid_scalar->concentration.col(i).data(),
                        n_species);
                    ++i;
                  });
  }

  void SimulationUnit::post_init_container(
      MC::DistributionVariantInt distribution_variant)
  {
    auto &to_process = mc_unit->container.to_process;

    const size_t n_compartments = mc_unit->domain.getNumberCompartments();

    //     #pragma omp target data map(to: to_process[0:to_process.size()])
    //     map(to: param) map(tofrom: mc_unit-&gt;domain)
    // {
    // #pragma omp target teams distribute parallel for schedule(static) \
//     map(to: n_compartments) \
//     num_threads(this-&gt;n_thread)
    //     for (auto it = to_process.begin(); it < to_process.end(); ++it)
    //     {
    //         auto prng = MC::PRNG::get_rng(omp_get_thread_num());
    //         auto distribution =
    //         MC::get_distribution_int&lt;size_t&gt;(param);

    //         auto &amp;&amp;particle = *it;
    //         particle.current_container = distribution(prng);

    //         auto &amp;i_container =
    //         mc_unit-&gt;domain[particle.current_container];

    //         kmodel.init_kernel(particle);
    //         __ATOM_INCR__(i_container.n_cells);
    //     }
    // }

#pragma omp parallel default(none),                                            \
    shared(to_process, n_compartments, distribution_variant)                   \
    num_threads(this->n_thread)
    {
      auto prng = MC::PRNG::get_rng(omp_get_thread_num());
      auto distribution =
          MC::get_distribution_int<size_t>(distribution_variant);
      const auto size_p = to_process.size();
#pragma omp for schedule(static)
      // for (auto it = to_process.begin(); it < to_process.end(); ++it)
      for (size_t i_p = 0; i_p < size_p; i_p++)
      {

        auto &&particle = to_process[i_p]; //*it;
        particle.current_container = distribution(prng);

        auto &i_container = mc_unit->domain[particle.current_container];

        kmodel.init_kernel(particle);
        __ATOM_INCR__(i_container.n_cells);
      }
    }
  }

  void SimulationUnit::cycleProcess(const double d_t)
  {

    const auto &cumulative_probability = flow_liquid->cumulative_probability;

    const auto &m_transition = this->flow_liquid->transition_matrix;
    const auto &_kmodel = this->kmodel;

    auto &unit = *this->mc_unit;
    auto contribs = this->liquid_scalar->getThreadContribs();
    auto &extras = this->mc_unit->extras;
    auto &to_process = this->mc_unit->container.to_process;
    const auto size = static_cast<size_t>(to_process.size());
#pragma omp for
    for (size_t i_particle = 0; i_particle < size; ++i_particle)
    {
      // particle->weight = 1./size;
      p_kernel(d_t,
               unit,
               contribs,
               extras,
               _kmodel,
               to_process[i_particle],
               m_transition,
               cumulative_probability);
    }

#pragma omp master
    post_process_reducing();
  }

  // #pragma omp for
  //     for (auto particle = to_process.begin(); particle < to_process.end();
  //          ++particle)
  //     {
  //       // particle->weight = 1./size;
  //       p_kernel(d_t,
  //                unit,
  //                contribs,
  //                extras,
  //                _kmodel,
  //                *particle,
  //                m_transition,
  //                cumulative_probability);
  //     }

  // #pragma omp master
  //     post_process_reducing();
  //   }

  void p_kernel(double d_t,
                MC::MonteCarloUnit &unit,
                std::span<Eigen::MatrixXd> _contribs,
                std::span<MC::ThreadPrivateData> _extras,
                const KModel &_kmodel,
                MC::Particles &particle,
                auto &m_transition,
                auto &cumulative_probability)
  {

    if (particle.status == MC::CellStatus::DEAD)
    {
      return;
    }

    auto &domain = unit.domain;

    const size_t i_thread = omp_get_thread_num();
    auto &events = unit.ts_events[i_thread];
    auto &thread_contrib = _contribs[i_thread];
    auto &thread_extra = _extras[i_thread];

    auto &rng = thread_extra.rng;
    const double rnd = rng.double_unfiform();
    const double rdn2 = rng.double_unfiform();

    const auto &concentrations =
        domain[particle.current_container].concentrations;

    kernel_move(
        rnd, rdn2, domain, particle, d_t, m_transition, cumulative_probability);

    _kmodel.update_kernel(d_t, particle, concentrations);

    if (particle.status == MC::CellStatus::DEAD)
    {
      events.incr<MC::EventType::Death>();
      __ATOM_DECR__(
          domain[particle.current_container].n_cells) // TODO: check overflow
      particle.clearState(MC::CellStatus::DEAD);
      thread_extra.in_dead_state.emplace_back(&particle);
    }
    else
    {
      if (particle.status == MC::CellStatus::CYTOKINESIS)
      {
        auto child = _kmodel.division_kernel(particle);
        events.incr<MC::EventType::NewParticle>();
        _kmodel.contribution_kernel(child, thread_contrib);
        __ATOM_INCR__(domain[child.current_container].n_cells)
        thread_extra.extra_process.emplace_back(std::move(child));
      }

      _kmodel.contribution_kernel(particle, thread_contrib);
    }
  };

} // namespace Simulation
