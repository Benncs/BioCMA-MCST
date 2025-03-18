#include <mc/unit.hpp>
#include <mc/mcinit.hpp>
auto init()
{

  const auto init = ExecInfo{
      .run_id = 1, .n_rank = 1, .current_rank = 0, .thread_per_process = 1, .verbose = false};
  const size_t n_particle = 10;
  std::vector<double> volumes = {20e-3};
  const double x0 = 1;
  double total_mass = 0;
  NeighborsView<HostSpace> neighb;
  {
    using current_model = DefaultModel;
    Kokkos::View<double**, Kokkos::LayoutRight, HostSpace> host_concentration(
        "test_host_concentration", 3, 1);
    host_concentration(0, 0) = 50.;
    auto device = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_concentration);
    LocalConcentrationView concentration = Kokkos::subview(device, Kokkos::ALL, 0);
    return MC::init<current_model>(n_particle, volumes, neighb, total_mass);
  }
}

int main()
{

  Kokkos::initialize();
  {
    auto unit = init();
  }


  Kokkos::finalize();
}
