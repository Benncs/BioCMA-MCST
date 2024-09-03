#ifndef __COMMON_KOKKOS_VECTOR__
#define __COMMON_KOKKOS_VECTOR__

#include <cmt_common/macro_constructor_assignment.hpp>
#include <Kokkos_Core.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <decl/Kokkos_Declare_OPENMP.hpp>
#include <impl/Kokkos_ViewCtor.hpp>
#include <utility>

using ComputeSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostSpace = Kokkos::HostSpace::memory_space;

template <typename T, typename Space> class KokkosVector
{
public:
  using Layout = Kokkos::LayoutRight;
  ~KokkosVector() = default;

  DEFAULT_COPY_MOVE_AC(KokkosVector)

  explicit KokkosVector(const size_t capacity,
                        bool alloc = true,
                        std::string label = "")
      : _owned_data(Kokkos::view_alloc(label + "_owned_data",
                                       Kokkos::WithoutInitializing))
  {
    __allocate__(capacity);
    if (alloc)
    {
      n_used_elements = capacity;
    }
    else
    {
      n_used_elements = 0;
    }
  }

  void clear()
  {
    n_used_elements = 0;
  }

  void reset()
  {
    n_used_elements = 0;
    n_allocated_element = 0;
    Kokkos::resize(_owned_data, n_allocated_element);
  }

  void resize(size_t n)
  {
    n_used_elements = 0;
    __allocate__(n);
  }

  auto static with_capacity(std::size_t capacity)
  {
    auto rhs = KokkosVector(capacity);
    rhs.n_used_elements = 0;
    return rhs;
  }

  KOKKOS_INLINE_FUNCTION auto &back()
  {
    return _owned_data(n_used_elements - 1);
  }
  KOKKOS_INLINE_FUNCTION auto &back() const
  {
    return _owned_data(n_used_elements - 1);
  }

  KOKKOS_INLINE_FUNCTION T &operator[](size_t i)
  {
    return _owned_data(i);
  }

  KOKKOS_INLINE_FUNCTION const T &operator[](size_t i) const
  {
    return _owned_data(i);
  }


  KOKKOS_INLINE_FUNCTION void set(size_t i,T&& data)const
  {
    _owned_data(i) = std::move(data);
  }

  template <typename MS1, typename MS2>
  static void migrate(KokkosVector<T, MS1> src, KokkosVector<T, MS2> &dest)
  {
    dest._owned_data = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                           src._owned_data);
    dest.n_used_elements = src.n_used_elements;
    dest.n_allocated_element = src.n_allocated_element;
  }

  KOKKOS_INLINE_FUNCTION auto size() const
  {
    return n_used_elements;
  }
  
  KOKKOS_INLINE_FUNCTION auto capacity() const
  {
    return n_allocated_element;
  }

  KOKKOS_FUNCTION bool emplace(T &&d)
  {
    const auto local_used_size = n_used_elements;
    if (local_used_size < n_allocated_element)
    {
      Kokkos::atomic_increment(&n_used_elements);
      _owned_data(local_used_size) = std::move(d);
      return true;
    }
    return false;
  }

  KOKKOS_INLINE_FUNCTION double get_allocation_factor()
  {
    return extra_allocation_factor;
  }

  KOKKOS_INLINE_FUNCTION void set_allocation_factor(double value)
  {
    assert(value >= 1);
    extra_allocation_factor = value;
  }

  void insert(const KokkosVector<T, Space> &rhs)
  {

    const auto original_size = n_used_elements;
    const auto n_add_item = rhs.size();
    __allocate__(original_size + rhs.size());

    auto data = this->data();
    Kokkos::parallel_for(
        "InsertNew",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_add_item),
        KOKKOS_LAMBDA(const size_t i) {
          data(original_size + i) = rhs._owned_data(i);
        });
    Kokkos::fence();

    n_used_elements += n_add_item;
  }

  Kokkos::View<T *, Layout, Space> _owned_data{};

protected:
  auto data()
  {
    return _owned_data;
  }
  auto ddata() const
  {
    return _owned_data;
  }

private:
  uint64_t n_allocated_element = 0;
  uint64_t n_used_elements = 0;
  double extra_allocation_factor = default_extra_allocation_factor;

  static constexpr double default_extra_allocation_factor = 1.5;

  auto __allocate__(const std::size_t new_size)
  {
    if (new_size > 0)
    {
      if (new_size >= n_allocated_element)
      {
        const auto new_allocated_size = static_cast<std::size_t>(
            std::ceil(static_cast<double>(new_size) * extra_allocation_factor));

        // std::cout << "New allocation from " << n_allocated_element << " to "
        //           << new_allocated_size << std::endl;

        n_allocated_element = new_allocated_size;
        Kokkos::resize(_owned_data, n_allocated_element);
      }
    }
    return new_size;
  }
};

#endif