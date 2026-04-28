#ifndef __CONTRIBUTION_KERNEL_HPP__
#define __CONTRIBUTION_KERNEL_HPP__
#include "Kokkos_Assert.hpp"
#include <Kokkos_Core.hpp>
#include <mc/particles_container.hpp>
#include <mc/traits.hpp>

template <ModelType M> struct ContributionFunctor
{
  struct Tag0D
  {
  };
  struct Tag3D
  {
  };

  using TeamPolicy = Kokkos::TeamPolicy<ComputeSpace>;
  using TeamMember = TeamPolicy::member_type;

  using ScratchView
      = Kokkos::View<float_t*,
                     TeamPolicy::execution_space::scratch_memory_space>;
  using float_t = float;

  ContributionFunctor(std::size_t particle_per_team,
                      MC::ContributionView contribution_scatter,
                      MC::ParticlesContainer<M> particles)
      : m_particle_per_team(particle_per_team),
        m_contribution_scatter(std::move(contribution_scatter)),
        m_particles(std::move(particles))
  {
  }

  size_t np{};

  void
  update(MC::ParticlesContainer<M> _particles)

  {
    this->m_particles = std::move(_particles);
    np = m_particles.n_particles();
  }

  std::size_t m_particle_per_team;

  MC::ContributionView m_contribution_scatter;
  MC::ParticlesContainer<M> m_particles;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(Tag0D _tag, const TeamMember& team) const
  {
    (void)_tag;
    const std::size_t p0 = team.league_rank() * m_particle_per_team;
    const std::size_t n_particle = m_particles.n_particles();
    const auto& status = m_particles.status;
    const auto& contribs = m_particles.contribs;
    static const auto n_c = M::n_c;

    ScratchView scratch(team.team_scratch(0), n_c);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, n_c),
                         [&](const std::size_t j)
                         { scratch(j) = float_t{ 0 }; });

    const auto upper_bound = ((p0 + m_particle_per_team) >= n_particle)
                                 ? n_particle - p0
                                 : m_particle_per_team;

    KOKKOS_ASSERT(upper_bound > 0 && upper_bound <= m_particle_per_team
                  && (p0 + upper_bound) <= n_particle);

    team.team_barrier();

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, 0, upper_bound),
        [&](const std::size_t relative_index)
        {
          const std::size_t flatten_index = p0 + relative_index;
          if (status(flatten_index) != MC::Status::Idle)
          {
            return;
          }
          const auto weight = m_particles.get_weight(flatten_index);
          for (std::size_t j = 0; j < n_c; ++j)
          {
            Kokkos::atomic_add(&scratch(j),
                               weight * contribs(flatten_index, j));
          }
        });

    team.team_barrier();
    const auto& cs = m_contribution_scatter;
    Kokkos::single(Kokkos::PerTeam(team),
                   [=]()
                   {
                     auto access = cs.access();
                     for (std::size_t j = 0; j < M::n_c; ++j)
                     {
                       access(j, 0) += scratch(j);
                     }
                   });
  }

  // KOKKOS_INLINE_FUNCTION
  // void
  // operator()(const Tag3D _tag, const TeamMember& team) const
  // {
  //   (void)_tag;

  //   const std::size_t p0 = team.league_rank() * m_particle_per_team;
  //   const auto& c = m_particles.contribs;
  //   const auto& status = m_particles.status;
  //   const auto n_particle = m_particles.n_particles();
  //   const auto& positions = m_particles.position;

  //   const auto upper_bound = ((p0 + m_particle_per_team) >= n_particle)
  //                                ? n_particle - p0 - 1
  //                                : m_particle_per_team;
  //   KOKKOS_ASSERT(upper_bound >= 0 && upper_bound < n_particle);

  //   constexpr std::size_t work_per_thread = 32;
  //   KOKKOS_ASSERT(m_particle_per_team%work_per_thread==0);

  //   Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, upper_bound),
  //                        [&](const std::size_t i)
  //                        {
  //                          const std::size_t p = p0 + i;
  //                          if (status(p) != MC::Status::Idle)
  //                          {
  //                            return;
  //                          }
  //                          auto access = m_contribution_scatter.access();
  //                          const double weight = m_particles.get_weight(p);
  //                          const auto pos = positions(p);

  //                          Kokkos::parallel_for(
  //                              Kokkos::ThreadVectorRange(team, 0, M::n_c),
  //                              [&](const int j)
  //                              { access(j, pos) += weight * c(p, j); });
  //                        });

  // }

  KOKKOS_INLINE_FUNCTION
  void
  operator()(const Tag3D _tag, const TeamMember& team) const
  {
    (void)_tag;

    const std::size_t p0 = team.league_rank() * m_particle_per_team;
    const auto& c = m_particles.contribs;
    const auto& status = m_particles.status;
    const auto n_particle = m_particles.n_particles();
    const auto& positions = m_particles.position;
    constexpr std::size_t work_per_thread = 32;

    const auto upper_bound = ((p0 + m_particle_per_team) >= n_particle)
                                 ? n_particle - p0
                                 : m_particle_per_team;
    KOKKOS_ASSERT(upper_bound >= 0 && upper_bound < n_particle);
    const auto range = (upper_bound + work_per_thread - 1) / work_per_thread;

    KOKKOS_ASSERT(m_particle_per_team % work_per_thread == 0);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, 0, range),
        [&](const std::size_t i)
        {
          auto access = m_contribution_scatter.access();
          for (std::size_t k = 0; k < work_per_thread; ++k)
          {
            const std::size_t p = p0 + i * work_per_thread + k;
            if (status(p) != MC::Status::Idle || i >= n_particle)
            {
              return;
            }
            const double weight = m_particles.get_weight(p);
            const auto pos = positions(p);
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, 0, M::n_c),
                                 [&](const int j)
                                 { access(j, pos) += weight * c(p, j); });
          }
        });
  }
};

#endif
