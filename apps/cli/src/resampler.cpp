

// TODO check case when serde is not enable
#include "Kokkos_Assert.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "impl/Kokkos_Profiling.hpp"
#include "mc/particles_container.hpp"
#include "mc/traits.hpp"
#include <core/case_data.hpp>
#include <iostream>
#include <mc/unit.hpp>
#include <stdexcept>
#include <variant>

template <ModelType Model>
void
copy_particle_to(MC::ParticlesContainer<Model> container1,
                 MC::ParticlesContainer<Model> container2,
                 size_t idx1,
                 size_t idx2)
{
  // KOKKOS_ASSERT(container.weights.extent(0) >= idx);
  KOKKOS_ASSERT(container1.status.extent(0) >= idx1);
  KOKKOS_ASSERT(container1.position.extent(0) >= idx1);

  KOKKOS_ASSERT(container2.status.extent(0) >= idx2);
  KOKKOS_ASSERT(container2.position.extent(0) >= idx2);

  // resampled_container.weights(i) = container.weights(idx);
  container1.position(idx1) = container2.position(idx2);
  container1.status(idx1) = container2.status(idx2);

  KOKKOS_ASSERT(container2.ages.extent(0) >= idx2);
  KOKKOS_ASSERT(container1.ages.extent(0) >= idx1);

  for (auto j = 0UL; j < container1.ages.extent(1); ++j)
  {
    container1.ages(idx1, j) = container2.ages(idx2, j);
  }

  KOKKOS_ASSERT(container2.contribs.extent(0) >= idx2);
  KOKKOS_ASSERT(container1.contribs.extent(0) >= idx1);

  for (auto j = decltype(Model::n_c){ 0 }; j < Model::n_c; ++j)
  {
    container1.contribs(idx1, j) = container2.contribs(idx2, j);
  }
  KOKKOS_ASSERT(container2.model.extent(0) >= idx2);
  KOKKOS_ASSERT(container1.model.extent(0) >= idx1);
  for (auto j = decltype(Model::n_var){ 0 }; j < Model::n_var; ++j)
  {
    container1.model(idx1, j) = container2.model(idx2, j);
  }
}

template <ModelType Model>
void
resampling(MC::ParticlesContainer<Model> container,
           std::size_t n_target,
           double initial_weight)
{
  if (n_target >= container.n_particles())
  {
    return;
  }

  MC::ParticlesContainer<Model> resampled_container(container.get_runtime(),
                                                    n_target);

  std::cout << "N particle  Before " << container.n_particles() << std::endl;
  std::cout << "N particle  After" << resampled_container.n_particles()
            << std::endl;
  auto n0 = container.n_particles();

  if constexpr (ConstWeightModelType<Model>)
  {

    double n1 = static_cast<double>(resampled_container.n_particles());
    double n2 = static_cast<double>(container.n_particles());

    // Unique weight then get only the first one
    double new_weight = n2 / n1 * container.get_weight(0);

    Kokkos::deep_copy(resampled_container.weights, new_weight);
  }
  else
  {
    throw std::runtime_error("Not implemented");
  }

  Kokkos::View<std::size_t*> kept_idx("kept_idx", n_target);

  Kokkos::Random_XorShift1024_Pool pool{ 1025 };

  Kokkos::parallel_for(
      "select_idx", Kokkos::RangePolicy(0, n_target), KOKKOS_LAMBDA(int i) {
        auto gen = pool.get_state();
        kept_idx(i) = gen.urand64(0, n0);
        pool.free_state(gen);
      });

  Kokkos::parallel_for(
      "swap_particles", Kokkos::RangePolicy(0, n_target), KOKKOS_LAMBDA(int i) {
        auto idx = kept_idx(i);
        copy_particle_to(resampled_container, container, i, idx);
      });

  // MC::ParticleWeigths<typename Model::FloatType> cumulative_weights(
  //     "cumulative_weight", container.capacity());

  // typename Model::FloatType total_weight{};
  // Kokkos::parallel_reduce(
  //     "k_total_weight",
  //     Kokkos::RangePolicy(0, container.n_particles()),
  //     KOKKOS_LAMBDA(int i, typename Model::FloatType& reducer) {
  //       reducer += container.weights(i) / initial_weight;
  //     },
  //     total_weight);

  // std::cout << total_weight << std::endl;
}

int
main()
{
  Kokkos::ScopeGuard g{};

  auto mc = Core::load_mc("/tmp/debug_serde_0.raw");

  if (mc == nullptr)
  {
    std::cout << "Cannot load" << std::endl;
    return 0;
  }

  const size_t n_target = 1000;
  const auto initial_weight = mc->init_weight;

  std::cout << "" << std::endl;

  auto callback = [n_target, initial_weight](auto&& container)
  { resampling(container, n_target, initial_weight); };

  std::visit(callback, mc->container);
}
