#include "Kokkos_Assert.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "decl/Kokkos_Declare_OPENMP.hpp"
#include <cassert>

#include <iostream>
#include <string>
#include <vector>

// A simple struct to test with KokkosVector
struct foo
{
  double a;
  int b;

  bool operator==(const foo& other) const
  {
    return a == other.a && b == other.b;
  }
};

void test_double()
{
  KokkosVector<double, Kokkos::HostSpace> vec(10);
  assert(vec.size() == 10);

  auto af = vec.get_allocation_factor();
  KOKKOS_ASSERT(af == 1.5);
  assert(vec.capacity() == vec.size() * af);

  vec.set_allocation_factor(2.);
  KOKKOS_ASSERT(vec.get_allocation_factor() == 2.);

  // Test adding elements
  vec.emplace(1.1);
  vec.emplace(2.2);
  vec.emplace(3.3);

  KOKKOS_ASSERT(vec.size() == 13);
  KOKKOS_ASSERT(vec[10] == 1.1);
  KOKKOS_ASSERT(vec[11] == 2.2);
  KOKKOS_ASSERT(vec[12] == 3.3);

  // Test back() and set()
  KOKKOS_ASSERT(vec.back() == 3.3);
  vec.set(1, 4.4);
  KOKKOS_ASSERT(vec[1] == 4.4);

  // Test clear and reset
  vec.clear();
  KOKKOS_ASSERT(vec.size() == 0);

  vec.reset();
  KOKKOS_ASSERT(vec.size() == 0);
  KOKKOS_ASSERT(vec.capacity() == 0);
}

void test_struct()
{
  // Create a KokkosVector of foo structs
  KokkosVector<foo, Kokkos::HostSpace> vec(10);

  // Test size and capacity
  assert(vec.size() == 10);
  assert(vec.capacity() == vec.size() * vec.get_allocation_factor());

  // Test adding elements
  vec.emplace({1.1, 2});
  vec.emplace({3.3, 4});

  assert(vec.size() == 12);
  assert(vec[10] == foo(1.1, 2));
  auto f = foo(3.3, 4);
  assert(vec[11] == f);

  // Test back() and set()
  assert(vec.back() == foo(3.3, 4));
  vec.set(1, foo(5.5, 6));
  assert(vec[1] == foo(5.5, 6));

  // Test clear and reset
  vec.clear();
  assert(vec.size() == 0);

  vec.reset();
  assert(vec.size() == 0);
  assert(vec.capacity() == 0);

  std::cout << "test_struct passed!\n";
}

void test_migration()
{
  // Create two KokkosVectors
  KokkosVector<double, Kokkos::HostSpace> src(10);
  KokkosVector<double, Kokkos::HostSpace> dest(10);

  src.emplace(1.1);
  src.emplace(2.2);

  // Migrate data from src to dest
  migrate(src, dest);

  // FIXME
  // assert(dest.size() == src.size());
  // assert(dest[0] == src[0]);
  // assert(dest[1] == src[1]);
}

void test_resize_and_capacity()
{
  // Create a KokkosVector
  KokkosVector<double, Kokkos::HostSpace> vec(5);

  // Resize the vector
  vec.resize(10);

  assert(vec.capacity() >= 10); // Capacity may be greater due to allocation factor
  assert(vec.size() == 0);

  std::cout << "test_resize_and_capacity passed!\n";
}

void in_kernel_emplace()
{
  KokkosVector<foo, Kokkos::DefaultExecutionSpace> vec =
      KokkosVector<foo, Kokkos::DefaultExecutionSpace>::with_capacity(15);

  auto lambda = KOKKOS_LAMBDA(const int i)
  {
    vec.emplace({5, i});
  };
    Kokkos::parallel_for(
        "test_in_kernel_emplace", vec.capacity(), lambda);



  auto l = Kokkos::create_mirror_view_and_copy(HostSpace(), vec._owned_data);

  
 
  std::cout << l.size() << std::endl;
  KOKKOS_ASSERT(l.size() >= 10); // Allocation factor 1.5

}

int main()
{
  Kokkos::initialize();
  {
    test_double();
    test_struct();
    test_migration();
    test_resize_and_capacity();
    in_kernel_emplace();
  }
  Kokkos::finalize();

  return 0;
}
