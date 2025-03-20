#include "Kokkos_Assert.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "biocma_cst_config.hpp"
#include "mc/traits.hpp"
#include <mc/particles_container.hpp>
#include <mc/unit.hpp>

template <ModelType M> void basic_test()
{
  const std::size_t size = 1000;
  MC::ParticlesContainer<M> container(size);

  KOKKOS_ASSERT(container.n_particles() == size)
  KOKKOS_ASSERT(container.capacity() == size * container.get_allocation_factor());
  KOKKOS_ASSERT(container.model.extent(0) == container.capacity());
  KOKKOS_ASSERT(container.model.extent(1) == M::n_var);

  KOKKOS_ASSERT(container.position.extent(0) == container.capacity());
  KOKKOS_ASSERT(container.status.extent(0) == container.capacity());

  if constexpr (ConstWeightModelType<M>)
  {
    KOKKOS_ASSERT(container.weights.extent(0) == 1);
  }
  else
  {
    KOKKOS_ASSERT(container.weights.extent(0) == container.capacity());
  }
}

template <ModelType M> void team_member_test()
{
  const std::size_t size = 1000;
  MC::ParticlesContainer<M> container(size);
  bool foo{};
  using TeamPolicy = Kokkos::TeamPolicy<ComputeSpace>;
  using TeamMember = TeamPolicy::member_type;
  std::size_t cumsum = 0;
  Kokkos::parallel_reduce(
      "spawn",
      MC::get_policty(foo, size),
      KOKKOS_LAMBDA(const TeamMember& team_handle, std::size_t& cs) {
        GET_INDEX(size);
        cs += 1;
      },
      cumsum);
  Kokkos::fence();
  KOKKOS_ASSERT(cumsum == size);
}

template <ModelType M> void div_test()
{
  const std::size_t size = 1000;
  const std::size_t ndiv = 10;
  MC::ParticlesContainer<M> container(size);
  MC::KPRNG::pool_type rng;
  Kokkos::parallel_for(
      "spawn", ndiv, KOKKOS_LAMBDA(const int i) {
        KOKKOS_ASSERT(container.handle_division(rng, i));
      });
  Kokkos::fence();
  KOKKOS_ASSERT(container.get_buffer_index() == ndiv);
}

template <ModelType M> void merge_test()
{
  const std::size_t size = 1000;
  const std::size_t ndiv = 10;
  MC::ParticlesContainer<M> container(size);
  MC::KPRNG::pool_type rng;
  Kokkos::parallel_for(
      "spawn", ndiv, KOKKOS_LAMBDA(const int i) {
        KOKKOS_ASSERT(container.handle_division(rng, i));
      });
  Kokkos::fence();
  KOKKOS_ASSERT(container.get_buffer_index() == ndiv);

  container.merge_buffer();

  KOKKOS_ASSERT(container.n_particles() == ndiv + size);
  KOKKOS_ASSERT(container.get_buffer_index() == 0);

  KOKKOS_ASSERT(container.model.extent(0) == container.position.extent(0) &&
                container.position.extent(0) == container.status.extent(0));

  if constexpr (ConstWeightModelType<M>)
  {
    KOKKOS_ASSERT(container.weights.extent(0) == 1);
  }
  else
  {
    KOKKOS_ASSERT(container.weights.extent(0) == container.status.extent(0));
  }
}

template <ModelType M> void clean_test()
{
  const std::size_t size = 1000;
  const std::size_t to_remove = 10;
  MC::ParticlesContainer<M> container(size);
  MC::KPRNG::pool_type rng(AutoGenerated::debug_MC_RAND_DEFAULT_SEED);
  Kokkos::parallel_for(
      "killparticle", to_remove, KOKKOS_LAMBDA(const int i) {
        (void)i;
        auto gen = rng.get_state();
        container.status(gen.urand64(0, to_remove)) = MC::Status::Dead;
        rng.free_state(gen);
      });
  Kokkos::fence();
    KOKKOS_ASSERT(container.n_particles() == size);
  container.clean_dead(to_remove);
  KOKKOS_ASSERT(container.n_particles() == size-to_remove);

}

template <ModelType M> void clean_test_and_shrink()
{
  const std::size_t size = 100;
  const std::size_t to_remove = 99;
  MC::ParticlesContainer<M> container(size);
  MC::KPRNG::pool_type rng(AutoGenerated::debug_MC_RAND_DEFAULT_SEED);
  Kokkos::parallel_for(
      "killparticle", to_remove, KOKKOS_LAMBDA(const int i) {
        (void)i;
        auto gen = rng.get_state();
        container.status(gen.urand64(0, to_remove)) = MC::Status::Dead;
        rng.free_state(gen);
      });
  Kokkos::fence();
    KOKKOS_ASSERT(container.n_particles() == size);
  container.clean_dead(to_remove);
  KOKKOS_ASSERT(container.n_particles() == size-to_remove);

}

int main()
{
  Kokkos::ScopeGuard guard;
  basic_test<DefaultModel>();
  team_member_test<DefaultModel>();
  div_test<DefaultModel>();
  merge_test<DefaultModel>();
  clean_test<DefaultModel>();
  clean_test_and_shrink<DefaultModel>();
}