
#ifndef __MC_MANIP_HPP__
#define __MC_MANIP_HPP__

#include <mc/particles_container.hpp>
#include <mc/traits.hpp>
namespace MC
{

  template <ModelType M>
  double
  total_mass(const MC::ParticlesContainer<M>& container)
  {

    double total_mass = 0.;
    Kokkos::parallel_reduce(
        "mc_total_mass",
        Kokkos::RangePolicy<>(0, container.n_particles()),
        KOKKOS_LAMBDA(const int i, double& red_mass) {
          const double w = container.get_weight(i);
          red_mass += M::mass(i, container.model) * w;
        },
        total_mass);
    KOKKOS_ASSERT(total_mass >= 0.);
    return total_mass;
  }

  template <ModelType M>
  double
  mean_mass(const MC::ParticlesContainer<M>& container)
  {

    double total_mass = 0.;
    double denum = 0.;
    Kokkos::parallel_reduce(
        "mc_total_mass",
        Kokkos::RangePolicy<>(0, container.n_particles()),
        KOKKOS_LAMBDA(const int i, double& red_mass, double& red_denum) {
          const double w = container.get_weight(i);
          red_denum += w;
          red_mass += M::mass(i, container.model) * w;
        },
        total_mass,
        denum);
    KOKKOS_ASSERT(total_mass >= 0.);
    KOKKOS_ASSERT(denum > 0.);
    return total_mass / denum;
  }

  template <ModelType M>
  double
  mean_age(const MC::ParticlesContainer<M>& container)
  {

    double total_age = 0.;
    double denum = 0.;
    Kokkos::parallel_reduce(
        "mc_total_mass",
        Kokkos::RangePolicy<>(0, container.n_particles()),
        KOKKOS_LAMBDA(const int i, double& red_a, double& red_denum) {
          const double w = container.get_weight(i);
          red_denum += w;
          red_a += w * container.ages(i, 1);
        },
        total_age,
        denum);
    KOKKOS_ASSERT(total_age >= 0.);
    KOKKOS_ASSERT(denum > 0.);
    return total_age / denum;
  }

} // namespace MC

#endif
