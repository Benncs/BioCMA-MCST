

// TODO check case when serde is not enable
#include "Kokkos_Assert.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "common/traits.hpp"
#include "impl/Kokkos_Profiling.hpp"
#include "mc/manip.hpp"
#include "mc/particles_container.hpp"
#include "mc/traits.hpp"
#include <core/case_data.hpp>
#include <iostream>
#include <mc/unit.hpp>
#include <memory>
#include <stdexcept>
#include <variant>

struct Args
{
  std::size_t n_target;
  double weight;
};

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

struct SamplingInvariants
{
  double total_mass;
  double mean_mass;
  double mean_age;
  std::vector<double> mean;

  template <ModelType M>
  explicit SamplingInvariants(const MC::ParticlesContainer<M>& container)
  {
    total_mass = MC::total_mass(container);
    mean_mass = MC::mean_mass(container);
    mean_age = MC::mean_age(container);
    mean = MC::mean_properties(container);
  }

  static bool
  validate(const SamplingInvariants& before, const SamplingInvariants& after)
  {

    std::cout << "Before" << std::endl;
    std::cout << "total mass: " << before.total_mass << std::endl;
    std::cout << "mean mass: " << before.mean_mass << std::endl;
    std::cout << "mean age: " << before.mean_age << std::endl;

    for (auto&& i : before.mean)
    {
      std::cout << "P :" << i << std::endl;
    }

    std::cout << "After " << std::endl;

    double var
        = (before.total_mass - after.total_mass) / (before.total_mass) * 100.;
    std::cout << "total mass: " << after.total_mass << "(" << var << "%)"
              << std::endl;
    var = (before.mean_mass - after.mean_mass) / (before.mean_mass) * 100.;
    std::cout << "mean mass: " << after.mean_mass << "(" << var << "%)"
              << std::endl;
    var = (before.mean_age - after.mean_age) / (before.mean_age) * 100.;
    std::cout << "mean age: " << after.mean_age << "(" << var << "%)"
              << std::endl;

    for (auto&& i : after.mean)
    {
      std::cout << "P :" << i << std::endl;
    }

    // FIXME: almost_equal is not enough, need to check rel and abs error
    // This way, rel can be fixed once for all invariants
    const double tol_1 = 1e-7;
    const double tol_2 = 1e-8;
    const double tol_3 = 1e-8;
    bool f = almost_equal(after.total_mass, before.total_mass, tol_1);
    f &= almost_equal(after.mean_mass, before.mean_mass, tol_2);
    f &= almost_equal(after.mean_age, before.mean_age, tol_3);

    return f;
  }

  bool
  operator==(const SamplingInvariants& rhs) const
  {

    bool f = rhs.total_mass == this->total_mass;

    f &= rhs.mean_mass == this->mean_mass;

    f &= rhs.mean_age == this->mean_age;

    return f;
  }
};

template <ModelType Model>
void
resampling(MC::ParticlesContainer<Model>& container, Args args)
{
  if (args.n_target >= container.n_particles())
  {
    return;
  }

  SamplingInvariants inv_1(container);

  MC::ParticlesContainer<Model> resampled_container(container.get_runtime(),
                                                    args.n_target);

  std::cout << "N particle before " << container.n_particles() << std::endl;
  std::cout << "N particle after " << resampled_container.n_particles()
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

  Kokkos::View<std::size_t*> kept_idx("kept_idx", args.n_target);

  Kokkos::Random_XorShift1024_Pool pool{ 1025 };

  Kokkos::parallel_for(
      "select_idx",
      Kokkos::RangePolicy(0, args.n_target),
      KOKKOS_LAMBDA(int i) {
        auto gen = pool.get_state();
        kept_idx(i) = gen.urand64(0, n0);
        pool.free_state(gen);
      });

  Kokkos::parallel_for(
      "swap_particles",
      Kokkos::RangePolicy(0, args.n_target),
      KOKKOS_LAMBDA(int i) {
        auto idx = kept_idx(i);
        copy_particle_to(resampled_container, container, i, idx);
      });

  SamplingInvariants inv_2(resampled_container);

  if (SamplingInvariants::validate(inv_1, inv_2))
  {
    std::cout << "Resample OK" << std::endl;
  }
  else
  {
    std::cout << "Resample does not conserve all invariants" << std::endl;
  }

  container = resampled_container;

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

  const auto serde_ref_path = "/tmp/debug_serde_0.raw";
  const auto new_serde_path = "/tmp/debug_rs_serde_0.raw";
  auto mc = Core::load_mc(serde_ref_path);

  std::cout << "" << std::endl;

  if (mc == nullptr)
  {
    std::cout << "Cannot load" << std::endl;
    return 0;
  }

  const size_t n_target = 100;
  const auto initial_weight = mc->init_weight;

  Args args{ n_target, initial_weight };

  auto callback = [args](auto& container) { resampling(container, args); };

  std::visit(callback, mc->container);

  Core::swap_mc(std::move(mc), serde_ref_path, new_serde_path);
}
