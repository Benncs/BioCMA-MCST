#include <kokkos_eigen.hpp>
using namespace KokkosEigen;

void
example_1()
{

  auto view_1 = KokkosEigen2D<float>("test", 6, 5);

  view_1.host_to_device_sync();

  auto view_2 = KokkosEigen2D<float,
                              Kokkos::LayoutRight,
                              Kokkos::DefaultExecutionSpace>("test2", 5, 4);
  view_2.host_to_device_sync();

  auto tmp = view_1.eigen() * view_2.eigen();
  tmp.eval();
  Kokkos::printf("%d %d\r\n", tmp.cols(), tmp.rows());
}
void
example_2()
{

  auto view_1 = KokkosEigen2D<float>("test", 5, 4);

  view_1.host_view()(0, 0) = 1;
  view_1.set_host(0, 1, 1.);
  assert(view_1.get_host(0, 1) == 1.);

  view_1.host_to_device_sync();

  auto view_2 = KokkosEigen2D<float, Kokkos::LayoutRight>("test2", 5, 4);
  view_2.host_to_device_sync();
  auto tmp = view_1.as_array() * view_2.as_array();

  Kokkos::printf("%d %d\r\n", tmp.cols(), tmp.rows());
}

int
main()
{
  Kokkos::ScopeGuard g{};
  {
    example_1();
    example_2();
  }
}
