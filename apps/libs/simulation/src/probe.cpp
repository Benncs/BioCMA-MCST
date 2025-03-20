// #include <Kokkos_Atomic.hpp>
// #include <Kokkos_Core.hpp>
// #include <Kokkos_StdAlgorithms.hpp>
// #include <common/kokkos_vector.hpp>
// #include <cstddef>
// #include <cstdint>
// #include <optional>
// #include <simulation/probe.hpp>


/* Old Implementation


#define __condition(__MACRO__)                                                 \
  if (this->n_p == 0 || this->n_t == 0)                                        \
  {                                                                            \
    __MACRO__                                                                  \
  }
#define TEST_INIT __condition(return;)
#define TEST_INIT_OPT __condition(return std::nullopt;)

namespace Simulation
{

  Probes::Probes(size_t initial_n_particle, size_t n_t_flush)
      : n_p(initial_n_particle), n_t(n_t_flush), n_s(0), labels("labels",
n_s),buffer("buffer"),internal_counter("i_c")
  {
    this->probes =
        Kokkos::View<double ***, ComputeSpace>("probes", n_p, n_t, n_s);
  }

  Probes::Probes(size_t initial_n_particle,
                 size_t n_t_flush,
                 std::span<std::string> _labels)
      : n_p(initial_n_particle), n_t(n_t_flush), n_s(_labels.size()),
        probes("probes", n_p, n_t, n_s), labels("labels",
n_s),buffer("buffer"),internal_counter("i_c")
  {

    for (size_t i = 0; i < n_s; ++i)
    {
      labels(i) = _labels[i];
    }
  }

  void Probes::set(double val) const
  {
    const uint64_t off = 1;
    auto i = Kokkos::atomic_fetch_add(&internal_counter(),off);
    if (i < buffer_size)
    {
      this->buffer(i) = val;
    }
  }

  double *Probes::raw_get()const
  {
    auto host = Kokkos::create_mirror_view_and_copy(HostSpace(),buffer);
    return host.data();
  }

  bool Probes::need_export() const
  {
    return internal_counter() >= buffer_size;
  }

  void Probes::clear()
  {
    Kokkos::deep_copy(buffer, 0);
    Kokkos::deep_copy(internal_counter, 0);
  }

  Probes::Probes(size_t initial_n_particle,
                 size_t n_t_flush,
                 std::initializer_list<std::string> _labels)
      : n_p(initial_n_particle), n_t(n_t_flush), n_s(_labels.size()),
        probes("probes", n_p, n_t, n_s), labels("labels", n_s),
buffer("buffer"),internal_counter("i_c")
  {

    size_t i = 0;
    for (auto &&s : _labels)
    {
      labels(i) = s;
    }
  }
  void Probes::add_probe(std::string_view label)
  {
    TEST_INIT
    n_s++;
    Kokkos::resize(probes, n_p, n_t, n_s);
    Kokkos::resize(labels, n_s);
    labels(n_s - 1) = label;
  }

  Probes::SubViewTimeType Probes::get()
  {
    return Kokkos::subview(probes, Kokkos::ALL(), get_counter++, Kokkos::ALL());
  }

  KOKKOS_FUNCTION std::optional<Probes::SubViewProbeType>
  Probes::get(std::string_view label)
  {

    TEST_INIT_OPT
    auto exespace = Kokkos::DefaultHostExecutionSpace();
    auto cbeg = Kokkos::Experimental::cbegin(labels);
    auto cend = Kokkos::Experimental::cend(labels);
    auto it = Kokkos::Experimental::find(exespace, cbeg, cend, label);

    if (it != cend)
    {
      size_t index = Kokkos::Experimental::distance(cbeg, it);
      return Kokkos::subview(probes, Kokkos::ALL(), Kokkos::ALL(), index);
    }
    return std::nullopt;
  }



  std::optional<std::array<size_t, 3>> Probes::get_raw(double **ptr)
  {

    TEST_INIT_OPT
    get_counter = 0;
    if (n_s != 0 && *ptr == nullptr)
    {
      host = Kokkos::create_mirror(HostSpace(), this->probes);
      Kokkos::deep_copy(host, this->probes);
      *ptr = host.data();

      return std::array<size_t, 3>({n_p, n_t, n_s});
    }
    return std::nullopt;
  }

  [[nodiscard]] std::span<std::string> Probes::get_labels() const
  {
    return {this->labels.data(), n_s};
  }
} // namespace Simulation

*/