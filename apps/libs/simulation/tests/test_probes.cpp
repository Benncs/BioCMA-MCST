#include <Kokkos_Core.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <simulation/probe.hpp>

void test_probes_set()
{
  auto probe = Simulation::Probes<2>(); // Buffer size 2
  assert(probe.set(1) == true);
  assert(probe.need_export() == false);
  assert(probe.set(2) == true);
  assert(probe.need_export() == true);
  assert(probe.set(3) == false);
}

void test_probes_get()
{
  constexpr std::size_t size_buff = 10;
  auto probe = Simulation::Probes<size_buff>(); // Buffer size 2
  for (auto i = 0LU; i < size_buff; ++i)
  {
    assert(probe.set(static_cast<double>(i)) == true);
  }
  assert(probe.set(55) == false);
  assert(probe.need_export() == true);
  auto rd = probe.get();
  assert(rd.size() == size_buff);

  for (auto i = 0LU; i < size_buff; ++i)
  {
    assert(rd[i] == static_cast<double>(i));
  }
}
void test_probes_clear()
{
  constexpr std::size_t size_buff = 20;

  constexpr double val_index_2 = 2. * 2.;

  auto probe = Simulation::Probes<size_buff>(); // Buffer size 2
  for (auto i = 0LU; i < size_buff; ++i)
  {
    assert(probe.set(2. * static_cast<double>(i)) == true); // NOLINT
  }
  assert(probe.need_export() == true);
  auto rd = probe.get();
  assert(rd.size() == size_buff);
  assert(rd[2] == val_index_2);
  probe.clear();
  assert(rd.size() == size_buff);
  assert(rd[2] != val_index_2);

  // Even if rd is deep copy of probe,svalue is reset to 0 but change doesn´t
  // change
  assert(probe.need_export() == false);
  for (auto i = 0LU; i < size_buff; ++i)
  {
    assert(probe.set(2. * static_cast<double>(i)) == true); // NOLINT
  }
  assert(probe.need_export() == true);
  rd = probe.get();
  assert(rd[2] == val_index_2);
}

int main()
{

  Kokkos::initialize();
  // test_probes_set();
  // std::cout<<"test_probes: Set OK"<<std::endl;
  // test_probes_get();
  // std::cout<<"test_probes: Get OK"<<std::endl;
  // test_probes_clear();
  // std::cout<<"test_probes: Clear OK"<<std::endl;
  Kokkos::finalize();
}

/*OLD IMPLEMENTATION TEST

#include <Kokkos_Assert.hpp>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <simulation/probe.hpp>
#include <span>
static constexpr uint64_t init_n_particle = 200;
static constexpr uint64_t n_time_flush = 10;

void test_default()
{
    Simulation::Probes p;
    p.add_probe("test");

}

void inner_test_basic(Simulation::Probes &p)
{
  auto opt = p.get("test_1");
  KOKKOS_ASSERT(opt.has_value());
  KOKKOS_ASSERT(opt->operator()(0, 0) == 50);
}


void test_basic()
{
  auto p = Simulation::Probes(init_n_particle, n_time_flush);
  p.add_probe("test_1");

  auto opt = p.get("test_1");

  KOKKOS_ASSERT(opt.has_value());
  auto probe = *opt;
  KOKKOS_ASSERT(probe.extent(0) == init_n_particle);
  KOKKOS_ASSERT(probe.extent(1) == n_time_flush);

  probe(0, 0) = 50;

  p.add_probe("test_2");
  inner_test_basic(p);
}

void test_basic_from_span()
{
  std::vector<std::string> labels = {"test_1", "test_2"};
  auto p = Simulation::Probes(init_n_particle, n_time_flush, labels);

  auto opt = p.get(labels[0]);

  KOKKOS_ASSERT(opt.has_value());
  auto probe = *opt;
  KOKKOS_ASSERT(probe.extent(0) == init_n_particle);
  KOKKOS_ASSERT(probe.extent(1) == n_time_flush);

  probe(0, 0) = 50;

  p.add_probe("test_3");
  inner_test_basic(p);
}

void test_ptr()
{
  {
    std::vector<std::string> labels = {"test_1", "test_2", "test_3"};
    auto probes = Simulation::Probes(init_n_particle, n_time_flush);
    for (auto &&i : labels)
    {
      probes.add_probe(i);
    }

    {
      auto probe = *probes.get(labels[0]); // Sure that it'll be valid
      probe(0, 0) = -5.;
    }
    {
      auto probe = *probes.get(labels[2]); // Sure that it'll be valid
      probe(0, 5) = M_PI;
    }

    auto res_label = probes.get_labels();
    KOKKOS_ASSERT(res_label.size() == labels.size());
    for (auto i = 0LU; i < res_label.size(); ++i)
    {
      KOKKOS_ASSERT(res_label[i] == labels[i]);
    }

    constexpr size_t n_sample = 3;
    double *ptr = nullptr;
    auto extents = probes.get_raw(&ptr);
    constexpr auto expected_dims =
        std::array<size_t, 3>({init_n_particle, n_time_flush, n_sample});

    KOKKOS_ASSERT(ptr != nullptr);
    KOKKOS_ASSERT(extents.has_value());
    KOKKOS_ASSERT(*extents == expected_dims);
    int i = 0; // First line
    int j = 0; // First row
    int k = 0; // First probe
    KOKKOS_ASSERT(ptr[i * (n_time_flush * n_sample) + j * n_sample + k] == -5.)
    i = 0;
    j = 5;
    k = 2;
    KOKKOS_ASSERT(ptr[i * (n_time_flush * n_sample) + j * n_sample + k] == M_PI)
  }

  {
    auto probes = Simulation::Probes(init_n_particle, n_time_flush);
    probes.add_probe("test_1");
    probes.add_probe("test_2");
    probes.add_probe("test_3");

    double a = 5;
    double *non_nullptr = &a;
    KOKKOS_ASSERT(!probes.get_raw(&non_nullptr).has_value());
  }

  {
    auto probes = Simulation::Probes(init_n_particle, n_time_flush);
    double *ptr = nullptr;
    KOKKOS_ASSERT(!probes.get_raw(&ptr).has_value())
  }
}

int main()
{

  Kokkos::initialize();
  test_default();
  test_basic();
  test_basic_from_span();
  test_ptr();
  Kokkos::finalize();
}

*/