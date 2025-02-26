
#include "common/kokkos_vector.hpp"
#include "embed.hpp"
#include "mc/particles/data_holder.hpp"
#include "mc/prng/prng.hpp"
#include "models/uptake.hpp"
#include <Kokkos_Core.hpp>
#include <core/post_process.hpp>
#include <mc/container_state.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_meta2.hpp>
#include <pybind11/embed.h>
#include <span>

using FuzzingModel = Models::Meta2;
using FuzzingParticle = MC::Particle<FuzzingModel>;

template <typename ExecSpace>
using ListType = Kokkos::View<FuzzingParticle**, Kokkos::LayoutRight, ExecSpace>;
using RngPool = Kokkos::Random_XorShift64_Pool<>;

static constexpr std::size_t n_fuzzing_particle = 10e4;
static constexpr size_t run_rng_seed = 1440;
static constexpr std::size_t n_fuzz_concentrations = 1;

template <typename F, Models::Uptake::UptakeModel<F> Model> void fuzz_uptake_init(Model& data)
{
  data.a_permease = 0.5;
  data.a_pts = 0.5;
  data.n_permease = 100;
}

template <typename GenType, typename F, Models::Uptake::UptakeModel<F> Model>
void fuzz_uptake_init2(GenType gen, Model& data)
{
  const double permease_dist = gen.drand(0., 1.); // Range for permease-related properties
  const double pts_dist = gen.drand(0., 1.);      // Range for permease-related properties
  const double n_permease_dist = gen.drand(0, 200); // Range for n_permease
  data.a_permease = permease_dist;
  data.a_pts = pts_dist;
  data.n_permease = n_permease_dist; // This one changes all
}

template <typename ExecSpace> ListType<ExecSpace> init(RngPool pool)
{
  ListType<ExecSpace> particles("particle", n_fuzzing_particle, n_fuzz_concentrations);
  MC::KPRNG globalrng;
  Kokkos::parallel_for(
      "Generate",
      Kokkos::RangePolicy<ExecSpace>(0, n_fuzzing_particle),
      KOKKOS_LAMBDA(const int i_particle) {
        auto gen = pool.get_state();

        {
          auto& p = particles(i_particle, 0);
          fuzz_uptake_init<float>(p.data);
          const double length_dist = gen.drand(1e-9, 5e-6); // Range for le
          p.data.length = length_dist;
          p.data.nu1 = gen.drand(0., 0.3) / (3600 * 1e3 * 1e12);
        }

        for (size_t i_fuzz = 1; i_fuzz < n_fuzz_concentrations; ++i_fuzz)
        {
          // auto& p2 = particles(i_particle, i_fuzz);
          // fuzz_uptake_init2<decltype(gen), float>(gen, p2.data);
          // const double length_dist = gen.drand(1e-9, 5e-6); // Range for le
          // p2.data.length = length_dist;
          // p2.data.nu1 = gen.drand(0., 0.3) / (3600 * 1e3 * 1e12);
          auto& p2 = particles(i_particle, i_fuzz);
          p2.init(globalrng);
        }
        pool.free_state(gen);
      });
  return particles;
}

template <typename ExecSpace> void fuzz(ListType<ExecSpace> particles, RngPool pool)
{
  const double d_t = 1e-3;
  const std::size_t n_compartment = 1;
  const std::size_t n_species = 4; // G,O,A,Co2
  Kokkos::View<double***, ExecSpace> environment(
      "concentrations", n_fuzz_concentrations, n_compartment, n_species);

  Kokkos::View<double***, ExecSpace> phi(
      "phi", n_fuzz_concentrations, n_species, n_fuzzing_particle);

  Kokkos::parallel_for(
      "fuzzing",
      Kokkos::RangePolicy<ExecSpace>(0, n_fuzz_concentrations),
      KOKKOS_LAMBDA(const int i_fuzz) {
        auto gen = pool.get_state();
        environment(i_fuzz, 0, 0) = gen.drand(0., 5.);
        environment(i_fuzz, 0, 1) = gen.drand(0., 9e-3);
        environment(i_fuzz, 0, 2) = gen.drand(0., 1);
        environment(i_fuzz, 0, 3) = 0;
        pool.free_state(gen);
      });
  Kokkos::fence();

  MC::KPRNG rng;
  for (std::size_t i_fuzz = 0; i_fuzz < n_fuzz_concentrations; ++i_fuzz)
  {
    auto subview = Kokkos::subview(environment, i_fuzz, 0, Kokkos::ALL);
    auto phi_subview = Kokkos::subview(phi, i_fuzz, Kokkos::ALL, Kokkos::ALL);

    Kokkos::parallel_for(
        "fuzzing",
        Kokkos::RangePolicy<ExecSpace>(0, n_fuzzing_particle),
        KOKKOS_LAMBDA(const int i_particle) {
          auto& p = particles(i_particle, i_fuzz);
          p.update(d_t, subview, rng);
          auto child = p.division(rng);
          phi_subview(0, i_particle) = p.data.nu1* (3600 * 1e3 * 1e12);
          phi_subview(1, i_particle) = child.data.nu1* (3600 * 1e3 * 1e12);
          
        });
    Kokkos::fence();

    // auto buffer = PostProcessing::get_properties(n_fuzzing_particle, particles, n_compartment);

    auto phi_host = Kokkos::create_mirror_view_and_copy(HostSpace(), phi_subview);
    auto c_host = Kokkos::create_mirror_view_and_copy(HostSpace(), subview);
    auto span_phi = std::span<double>(phi_host.data(), phi_host.size());
    auto span_c = std::span<double>(c_host.data(), c_host.size());

    EmbedPython::call_add_function(span_phi, n_species, n_fuzzing_particle, span_c, 1, n_species);
  }

  EmbedPython::call_show_function();
}

int main()
{
  Kokkos::initialize();
  pybind11::scoped_interpreter guard{};

  using Space = Kokkos::DefaultExecutionSpace;
  {
    RngPool pool(run_rng_seed);
    auto particles = init<Space>(pool);
    fuzz<Space>(particles, pool);
  }

  Kokkos::finalize();
}