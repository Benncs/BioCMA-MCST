
#include <core/post_process.hpp>
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <iostream>
#include <mc/container_state.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_ecoli.hpp>

using FuzzingModel = Models::Ecoli;
using FuzzingParticle = MC::Particle<FuzzingModel>;

template <typename ExecSpace>
using ListType = Kokkos::View<FuzzingParticle*, Kokkos::LayoutRight, ExecSpace>;
using RngPool = Kokkos::Random_XorShift64_Pool<>;

static constexpr std::size_t n_fuzzing_particle = 100;
static constexpr size_t run_rng_seed = 1440;

MC::ContainerState get_container(LocalConcentrationView concentrations)
{
  MC::ContainerState container;
  container.concentrations = concentrations;
  container.id = 0;
  container.n_cells = n_fuzzing_particle;
  container.volume_liq = 20e-3;
  container.volume_gas = 0.;
  return container;
}

template <typename ExecSpace> ListType<ExecSpace> init(RngPool pool)
{
  ListType<ExecSpace> particles("particle", n_fuzzing_particle);

  Kokkos::parallel_for(
      "Generate",
      Kokkos::RangePolicy<ExecSpace>(0, n_fuzzing_particle),
      KOKKOS_LAMBDA(const int i_particle) {
        auto gen = pool.get_state();
        const double permease_dist = gen.drand(0., 1.);   // Range for permease-related properties
        const double pts_dist = gen.drand(0., 1.);        // Range for permease-related properties
        const double length_dist = gen.drand(1e-9, 5e-6); // Range for le
        const double n_permease_dist = gen.drand(0, 200); // Range for n_permease
        pool.free_state(gen);
        auto& p = particles(i_particle);
        p.data.a_permease = permease_dist;
        p.data.a_pts = pts_dist;
        p.data.lenght = length_dist;
        p.data.n_permease = n_permease_dist;
        ;
      });
  return particles;
}

template <typename ExecSpace> void fuzz(ListType<ExecSpace> particles, RngPool pool)
{
  const double d_t = 1e-3;

  const std::size_t n_fuzz_concentrations = 50;
  const std::size_t n_compartment = 1;
  const std::size_t n_species = 4; // G,O,A,Co2
  Kokkos::View<double***, ExecSpace> environment(
      "concentrations", n_compartment, n_fuzz_concentrations, n_species);

  Kokkos::parallel_for(
      "fuzzing",
      Kokkos::RangePolicy<ExecSpace>(0, n_fuzz_concentrations),
      KOKKOS_LAMBDA(const int i_fuzz) {
        environment(0, i_fuzz, 0) = 0;
        environment(0, i_fuzz, 1) = 0;
        environment(0, i_fuzz, 2) = 0;
        environment(0, i_fuzz, 3) = 0;
      });
  Kokkos::fence();

  for (std::size_t i_fuzz = 0; i_fuzz < n_fuzz_concentrations; ++i_fuzz)
  {
    std::cout << "Fuzzing concentration: " << i_fuzz << std::endl;
    auto subview = Kokkos::subview(environment, 1, i_fuzz, Kokkos::ALL);
    Kokkos::parallel_for(
        "fuzzing",
        Kokkos::RangePolicy<ExecSpace>(0, n_fuzzing_particle),
        KOKKOS_LAMBDA(const int i_particle) {
          Kokkos::printf("Fuzzing particle %ld\r\n", i_particle);
          auto& p = particles(i_particle);
          p.update(d_t, subview, pool);
        });
    Kokkos::fence();

    auto buffer = PostProcessing::get_properties(n_fuzz_concentrations, particles, n_compartment);

    std::cout<<buffer.particle_values(0,0)<<std::endl;
  }
}

int main()
{
  Kokkos::initialize();

  using Space = Kokkos::DefaultExecutionSpace;
  {
    RngPool pool(run_rng_seed);
    auto particles = init<Space>(pool);
    fuzz<Space>(particles, pool);
  }
  Kokkos::finalize();
}