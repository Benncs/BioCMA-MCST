#include "Kokkos_Macros.hpp"
#include "impl/Kokkos_Profiling.hpp"
#include <Kokkos_Assert.hpp>
#include <kokkos_eigen.hpp>

using namespace KokkosEigen;

template <typename layout>
void
_impl_test()
{
  const std::size_t n_row = 10;
  const std::size_t n_col = 6;
  auto ke = KokkosEigen2D<float, layout>("test", n_row, n_col);

  Kokkos::printf("%d %d\r\n", ke.n_col(), ke.n_row());
  KOKKOS_ASSERT(ke.n_col() == n_col);
  KOKKOS_ASSERT(ke.n_row() == n_row);

  auto& eigen = ke.eigen();
  KOKKOS_ASSERT(eigen.cols() == n_col);
  KOKKOS_ASSERT(eigen.rows() == n_row);

  eigen.coeffRef(2, 3) = 4.;
  KOKKOS_ASSERT(ke.get_host(2, 3) == 4.);

  ke.set_host(2, 3, 5.);
  KOKKOS_ASSERT(eigen.coeff(2, 3) == 5.);
}

int
main()
{
  Kokkos::ScopeGuard g{};

  _impl_test<Kokkos::LayoutLeft>();
  _impl_test<Kokkos::LayoutRight>();
}
