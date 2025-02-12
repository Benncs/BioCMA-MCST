
#include "Kokkos_Core.hpp"
#include "common/kokkos_vector.hpp"
#include "decl/Kokkos_Declare_OPENMP.hpp"
#include "mc/mcinit.hpp"
#include "mc/particles/particle_list.hpp"
#include "mc/prng/prng.hpp"
#include "models/model_monod.hpp"
#include <Kokkos_Assert.hpp>
#include <mc/unit.hpp>

#include <mc/particles/particle_model.hpp>

template <typename ListType, typename ExtraType> struct TestKernel
{

  TestKernel(ListType _list, ExtraType _extra, auto _concentraiton, auto _rng)
      : list(std::move(_list)), extra(std::move(_extra)), concentration(std::move(_concentraiton)),
        rng(std::move(_rng))
  {
  }

  KOKKOS_FUNCTION void operator()(const std::size_t i_particle) const
  {
    auto& particle = list._owned_data(i_particle);
    KOKKOS_ASSERT(particle.properties.id == i_particle);
    particle.update(0.1, concentration, rng);
    const auto* np = extra.spawn();
    KOKKOS_ASSERT(np != nullptr);
    Kokkos::printf("%d\r\n", extra.size());
  }

  ListType list;
  ExtraType extra;
  LocalConcentrationView concentration;
  MC::KPRNG rng;
};

int main()
{
  Kokkos::initialize();
  const auto init = ExecInfo{1, 1, 0, 1, false};
  const size_t n_particle = 10;
  std::vector<double> volumes = {20e-3};
  const double x0 = 1;
  double total_mass = 0;
  CmaRead::Neighbors::Neighbors_const_view_t neighb;
  {
    using current_model = Models::Monod;
    Kokkos::View<double**, Kokkos::LayoutRight, HostSpace> host_concentration(
        "test_host_concentration", 3, 1);
    host_concentration(0, 0) = 50.;
    auto device = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_concentration);
    LocalConcentrationView concentration = Kokkos::subview(device, Kokkos::ALL, 0);
    auto unit = MC::init<current_model>(init, n_particle, volumes, neighb, x0, total_mass);
    auto rng = unit->rng;
    auto& container = std::get<MC::ParticlesContainer<current_model>>(unit->container);
    auto list = container.get_compute();

    MC::ParticleList<ComputeSpace, current_model> extra = container.process_buffer;

    Kokkos::printf("init %d\r\n", extra.size());

    for (int i = 0; i < 3; ++i)
    {
      auto k = TestKernel(list, extra, concentration, rng);
      Kokkos::parallel_for(n_particle, k);
    }
    Kokkos::fence();
    Kokkos::printf("final ??%d\r\n", extra.size());
    unit.reset();
    std::cout << "OK" << std::endl;
  }
  Kokkos::finalize();
  return 0;
}
