#include "Kokkos_Assert.hpp"
#include "Kokkos_Core_fwd.hpp"
#include <Kokkos_Core.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <iostream>
#include <cereal/details/traits.hpp>
#include <common/has_serialize.hpp>


void test_1d()
{
  const std::size_t n = 10;
  Kokkos::View<int*> view("view", n);
  Kokkos::parallel_for("initialize", n, KOKKOS_LAMBDA(const int i) { view(i) = i; });
  Kokkos::fence();

  std::ostringstream buff(std::ios::binary);
  cereal::BinaryOutputArchive oarchive(buff);
  serialize_view(oarchive, view);
  std::istringstream iss(buff.str(), std::ios::binary);
  cereal::BinaryInputArchive iarchive(iss);

  Kokkos::View<int*> new_view;
  deserialize_view(iarchive, new_view);
  assert(new_view.extent(0) == n && "size");
  Kokkos::parallel_for("verify",
                       n,
                       KOKKOS_LAMBDA(const int i){KOKKOS_ASSERT(new_view(i) == i)

                       });
  Kokkos::fence();
}

void test_2d()
{
  constexpr std::size_t n_p = 4;
  const std::size_t n = 10;
  Kokkos::View<float**> view("view", n, n_p);
  Kokkos::parallel_for(
      "initialize", n, KOKKOS_LAMBDA(const int i) {
        for (int j = 0; j < n_p; ++j)
        {
          view(i, j) = float(i) * float(j);
        }
      });
  Kokkos::fence();

  std::ostringstream buff(std::ios::binary);
  cereal::BinaryOutputArchive oarchive(buff);
  serialize_view(oarchive, view);
  std::istringstream iss(buff.str(), std::ios::binary);
  cereal::BinaryInputArchive iarchive(iss);
  decltype(view) new_view;
  deserialize_view(iarchive, new_view);

  assert(new_view.extent(0) == n && "row");
  assert(new_view.extent(1) == n_p && "col");
  Kokkos::parallel_for(
      "verify", n, KOKKOS_LAMBDA(const int i) {
        for (int j = 0; j < n_p; ++j)
        {
          KOKKOS_ASSERT(Kokkos::abs(new_view(i, j) - float(i) * float(j)) < 1e-6);
        }
      });
  Kokkos::fence();
}

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    test_1d();
    test_2d();
  }
  Kokkos::finalize();
  return 0;
}

