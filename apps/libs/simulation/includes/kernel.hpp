#ifndef __SIMULATION_KERNEL_HPP__
#define __SIMULATION_KERNEL_HPP__

#include "cma_read/light_2d_view.hpp"
#include "common/common.hpp"
#include "common/thread_safe_container.hpp"
#include "mc/events.hpp"
#include "mc/thread_private_data.hpp"
#include "models/types.hpp"
#include "transport.hpp"
#include <Eigen/Core>
#include <Kokkos_Core.hpp>
#include <Kokkos_Printf.hpp>
#include <Kokkos_Random.hpp>
#include <cstdio>
#include <mc/container_state.hpp>
#include <mc/particles/mcparticles.hpp>

struct Kernel
{

  using value_type = MC::Results;

  std::span<MC::Particles> to_process;
  KModel kmodel;
  Kokkos::Random_XorShift1024_Pool<> random_pool;
  double d_t{};
  std::span<const double> diag;
  CmaRead::L2DView<const double> view_cumulative_probability;
  Kokkos::View<MC::ContainerState *, Kokkos::LayoutRight> domain;
  Kokkos::View<const size_t **, Kokkos::LayoutStride> neighbors_view;
  Kokkos::View<double **, Kokkos::LayoutLeft> contribs;

  // Kokkos::View<MC::Particles*, Kokkos::LayoutLeft> particles;

  Kernel(auto model, auto random_pool, auto _domain_view)
      : kmodel(model), random_pool(random_pool),
        domain(std::move(_domain_view)){};

  ~Kernel() = default;
  Kernel() = default;
  Kernel(const Kernel &) = default;
  Kernel(Kernel &&) noexcept = default;
  Kernel &operator=(const Kernel &) = default;
  Kernel &operator=(Kernel &&) noexcept = default;

  void update(auto _to_process,
              double _d_t,
              auto _diag,
              auto _view_cumulative_probability,
              auto _neighbors_view,
              auto _contribs)
  {
    to_process = _to_process;
    d_t = _d_t;
    diag = _diag;
    view_cumulative_probability = _view_cumulative_probability;
    neighbors_view = _neighbors_view;
    contribs = _contribs;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(auto &&i_particle, MC::Results &results) const
  {

    auto &particle = to_process[i_particle];
    auto &status = particle.status;
    if (status == MC::CellStatus::DEAD)
    {
      return;
    }

    Eigen::MatrixXd local_contribs(contribs.extent(0), 1); // FIXME
    local_contribs.setZero();

    auto generator = random_pool.get_state();
    const double random_number_1 = generator.drand(0., 1.);
    const double random_number_2 = generator.drand(0., 1.);
    const double random_number_3 = generator.drand(0., 1.);
    random_pool.free_state(generator);

    Simulation::kernel_move(random_number_1,
                            random_number_2,
                            particle,
                            d_t,
                            diag,
                            view_cumulative_probability,
                            domain,
                            neighbors_view);

    Simulation::kernel_exit(d_t, random_number_3, domain, particle);

    if (status == MC::CellStatus::OUT)
    {
      results.events.incr<MC::EventType::Exit>();
      __ATOM_DECR__(domain[particle.current_container].n_cells)
      particle.clearState(MC::CellStatus::DEAD);
      results.index_in_dead_state.emplace_back(i_particle);
      return;
    }

    kmodel.update_kernel(
        d_t, particle, domain[particle.current_container].concentrations);

    if (status == MC::CellStatus::DEAD)
    {
      results.events.incr<MC::EventType::Death>();
      __ATOM_DECR__(domain[particle.current_container].n_cells)
      particle.clearState(MC::CellStatus::DEAD);
      return;
    }

    if (status == MC::CellStatus::CYTOKINESIS)
    {
      results.events.incr<MC::EventType::NewParticle>();
      status = MC::CellStatus::IDLE;

      __ATOM_INCR__(domain[particle.current_container].n_cells);
      auto child = kmodel.division_kernel(particle);
      kmodel.contribution_kernel(child, local_contribs);

      results.extra_process.emplace_back(std::move(child));
    }

    kmodel.contribution_kernel(particle, local_contribs);

    for (size_t i_species = 0; i_species < contribs.extent(0); ++i_species)
    {
      Kokkos::atomic_add(&contribs(i_species, particle.current_container),
                         local_contribs(EIGEN_INDEX(i_species), 0));
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(MC::Results &dst, const MC::Results &src) const
  {
    dst.events.inplace_reduce(src.events);

    dst.extra_process.insert(dst.extra_process.end(),
                             src.extra_process.begin(),
                             src.extra_process.end());
    dst.index_in_dead_state.insert(dst.index_in_dead_state.end(),
                                   src.index_in_dead_state.begin(),
                                   src.index_in_dead_state.end());
  }

  // KOKKOS_INLINE_FUNCTION
  // void init(MC::Results &val) const
  // {
  //   MC::init_r(val);
  // }
};

#endif
