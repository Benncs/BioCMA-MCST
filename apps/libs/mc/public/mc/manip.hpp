
#ifndef __MC_MANIP_HPP__
#define __MC_MANIP_HPP__

#include "Kokkos_Core.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "common/common.hpp"
#include "impl/Kokkos_CheckUsage.hpp"
#include "impl/Kokkos_Profiling.hpp"
#include "mc/alias.hpp"
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

  // TODO: instead of mean, calculate more moments   or quantile
  template <ModelType M>
  std::vector<double>
  mean_properties(const MC::ParticlesContainer<M>& container)
  {

    Kokkos::View<double[M::n_var]> means("means_properties");
    std::vector<double> vec_means;
    vec_means.reserve(M::n_var);

    double denum = 0.;
    const double overflow_limit = 1e30;
    auto model = container.model;

    auto policy = Kokkos::MDRangePolicy<ComputeSpace, Kokkos::Rank<2> >(
        { 0, 0 }, { container.n_particles(), M::n_var });

    Kokkos::parallel_reduce(
        "total_weights",
        Kokkos::RangePolicy<>(0, container.n_particles()),
        KOKKOS_LAMBDA(const int i, double& red) {
          red += container.get_weight(i) / overflow_limit;
        },
        denum);

    Kokkos::parallel_for(
        "mc_means", policy, KOKKOS_LAMBDA(const int i, const int j) {
          const double w = container.get_weight(i);
          const auto v = model(i, j);
          Kokkos::atomic_add(&means(j), v * w);
        });

    const auto tmp = Kokkos::create_mirror_view_and_copy(HostSpace(), means);
    for (int i = 0; i < M::n_var; ++i)
    {
      double real_mean = tmp(i) / (overflow_limit * denum);
      vec_means.push_back(real_mean);
    }

    return vec_means;
  }

} // namespace MC

#endif
