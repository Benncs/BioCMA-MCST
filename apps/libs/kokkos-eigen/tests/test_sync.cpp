#include <Kokkos_Assert.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <kokkos_eigen.hpp>

using namespace KokkosEigen;

template <typename layout>
void
_impl_test()
{
  const std::size_t n_row = 10;
  const std::size_t n_col = 6;
  auto ke = KokkosEigen2D<float, layout>("test", n_row, n_col);

  auto& eigen = ke.eigen();

  auto dev = ke.device_view();
  Kokkos::parallel_for(
      "_test_device_host", n_row, KOKKOS_LAMBDA(const int i_row) {
        dev(i_row, 0) = (float)i_row;
      });
  Kokkos::fence();
  ke.device_to_host_sync();

  for (std::size_t i = 0; i < n_row; ++i)
  {
    Kokkos::printf("Got %f, expected %f \r\n", eigen.coeff(i, 0), (float)i);
    KOKKOS_ASSERT(eigen.coeff(i, 0) == (float)i);
    KOKKOS_ASSERT(ke.get_host(i, 0) == (float)i);
  }

  for (std::size_t i = 0; i < n_col; ++i)
  {
    eigen.coeffRef(0, i) = i * 2;
  }
  ke.host_to_device_sync();
  Kokkos::parallel_for(
      "_test_host_device", n_col, KOKKOS_LAMBDA(const int i_col) {
        KOKKOS_ASSERT(dev(0, i_col) == 2 * i_col);
      });
}

int
main()
{
  Kokkos::ScopeGuard g{};

  _impl_test<Kokkos::LayoutLeft>();
  _impl_test<Kokkos::LayoutRight>();
}
