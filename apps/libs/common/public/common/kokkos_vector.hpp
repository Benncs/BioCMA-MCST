#ifndef __COMMON_KOKKOS_VECTOR__
#define __COMMON_KOKKOS_VECTOR__

#include "Kokkos_Core_fwd.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

using ComputeSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostSpace = Kokkos::HostSpace::memory_space;

/**
 * @brief A basic wrapper around Kokkos View to handle non-fixed size containers
 * with an API compatible with `std::vector`.
 *
 * This class provides a dynamic array-like container built on top of Kokkos
 * View, allowing it to manage memory in a parallel-friendly way while offering
 * an interface similar to the standard C++ `std::vector`. It supports
 * operations like resizing, inserting elements, and migrating data between
 * different memory spaces.
 *
 * Like all other dynamic container, the real amount of allocated memory is
 * bigger than one requested by user to avoid reallocating memory when a single
 * element is pushed into vector
 *
 * As this is a wrapper around kokkos_view, it's a non owning type therefore,
 * copying this class performs a shallow copy with an atomic increment
 *
 * @tparam T The type of elements stored in the vector.
 * @tparam Space The memory space where the data is stored (e.g., host or device
 * memory).
 */
template <typename T, typename Space> class KokkosVector
{
public:
  /// The layout used for the underlying Kokkos View (LayoutRight by default)
  using Layout = Kokkos::LayoutRight;
  /// The Kokkos View holding the data elements.
  Kokkos::View<T*, Layout, Space, Kokkos::MemoryTraits<Kokkos::Restrict>> _owned_data{};
  /**
   * @brief Destructor for KokkosVector.
   */
  ~KokkosVector() = default;

  /**
   * @brief Default copy and move constructors and assignment operators.
   */

  KokkosVector() = default;
  KokkosVector(const KokkosVector&) = default;
  KokkosVector(KokkosVector&&) = default;
  KokkosVector& operator=(const KokkosVector&) = default;
  KokkosVector& operator=(KokkosVector&&) = default;

  /**
   * @brief Constructs a KokkosVector with a specified capacity.
   *
   * @param capacity The initial capacity of the vector.
   * @param alloc If true, set the size as full if not memory is already
   * allocated push vector has size 0
   * @param label A string label for the Kokkos View used internally.
   */
  explicit KokkosVector(size_t capacity, bool alloc = true, std::string label = "");

  /**
   * @brief Clears the vector, setting the used elements count to zero.
   *  @warning Cannot be called inside kernel
   */
  void clear() noexcept;

  /**
   * @brief Resets the vector, clearing all elements and deallocating memory.
   *  @warning Cannot be called inside kernel
   */
  void reset();

  /**
   * @brief Resizes the vector to the specified size.
   *
   * @warning Cannot be called inside kernel
   * @param n The new size of the vector.
   */
  void resize(size_t n);

  /**
   * @brief Creates a KokkosVector with the specified capacity, without
   * allocating elements immediately.
   *
   * @warning Cannot be called inside kernel
   * @param capacity The capacity of the vector.
   * @return A KokkosVector instance with the specified capacity.
   */
  static KokkosVector<T, Space> with_capacity(std::size_t capacity);

  /**
   * @brief Accesses the last element in the vector.
   *
   * @return Reference to the last element.
   */
  KOKKOS_INLINE_FUNCTION auto& back()
  {
    return _owned_data(n_used_elements() - 1);
  }
  /**
   * @brief Accesses the last element in the vector (const version).
   *
   * @return Const reference to the last element.
   */
  KOKKOS_INLINE_FUNCTION auto& back() const
  {
    return _owned_data(n_used_elements() - 1);
  }

  /**
   * @brief Accesses an element by index.
   *
   * @param i The index of the element.
   * @return Reference to the element at index `i`.
   */
  KOKKOS_INLINE_FUNCTION T& operator[](size_t i)
  {
    KOKKOS_ASSERT(i < n_used_elements()); // Kokkos view is not bound checked when use release
    return _owned_data(i);
  }

  /**
   * @brief Accesses an element by index (const version).
   *
   * @param i The index of the element.
   * @return Const reference to the element at index `i`.
   */
  KOKKOS_INLINE_FUNCTION const T& operator[](size_t i) const
  {
    // Kokkos view is not bound checked when use release
    KOKKOS_ASSERT(i < n_used_elements());
    return _owned_data(i);
  }

  /**
   * @brief Sets an element at a specific index.
   *
   * This function is marked const because underlying data has internal
   * mutability
   * @param i The index at which to set the element.
   * @param data The element to set.
   */
  KOKKOS_INLINE_FUNCTION void set(size_t i, T&& data) const
  {
    // Kokkos view is not bound checked when use release
    KOKKOS_ASSERT(i < n_used_elements());

    _owned_data(i) = std::move(data);
  }

  /**
   * @brief Migrates the data from one KokkosVector to another.
   *
   * @warning Cannot be called inside kernel
   * @tparam MS1 The memory space of the source vector.
   * @tparam MS2 The memory space of the destination vector.
   * @param src The source vector.
   * @param dest The destination vector.
   */
  template <typename MS1, typename MS2>
  static void migrate(KokkosVector<T, MS1> src, KokkosVector<T, MS2>& dest)
  {
    dest._owned_data = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), src._owned_data);
    //TODO 
    // auto local_n_used_elements = src.n_used_elements;
    // auto dest_used = Kokkos::create_mirror_view_and_copy(Kokkos::SharedSpace(),local_n_used_elements);
    // dest.n_used_elements = dest_used;
    // dest.n_allocated_element = src.n_allocated_element;
  }

  /**
   * @brief Gets the number of elements currently in use.
   *
   * @return The number of used elements.
   */
  KOKKOS_INLINE_FUNCTION auto size() const
  {
    return n_used_elements();
  }

  /**
   * @brief Gets the capacity of the vector.
   *
   * @return The total capacity of the vector.
   */
  KOKKOS_INLINE_FUNCTION auto capacity() const
  {
    return n_allocated_element;
  }

  /**
   * @brief Attempts to add a new element to the vector.
   *
   * @param d The element to add.
   * @return True if the element was added successfully, false if the vector is
   * full.
   */
  KOKKOS_FUNCTION bool emplace(T&& d) const;

  /**
   * @brief Gets the current allocation factor.
   *
   * @warning Cannot be called inside kernel
   * @return The allocation factor used for resizing the vector.
   */
  [[nodiscard]] double get_allocation_factor() const noexcept;

  /**
   * @brief Sets the allocation factor.
   *
   * @warning Cannot be called inside kernel
   * @param value The new allocation factor (must be >= 1).
   */
  void set_allocation_factor(double value) noexcept;

  /**
   * @brief Inserts the elements from another KokkosVector into this vector.
   *
   * @warning Cannot be called inside kernel
   * @param rhs The vector whose elements will be inserted.
   */
  void insert(const KokkosVector<T, Space>& rhs);

  /**
   * @brief Gets the underlying data.
   *
   * @return A Kokkos View of the data.
   */
  KOKKOS_INLINE_FUNCTION auto data()
  {
    return _owned_data;
  }

  template <class Archive> void save(Archive& ar) const
  {
    auto host_particle = Kokkos::create_mirror_view_and_copy(HostSpace(), this->_owned_data);
    std::vector<T> data_vector(host_particle.data(), host_particle.data() + n_used_elements());
    ar(n_allocated_element, n_used_elements(), extra_allocation_factor, data_vector);
  }

  template <class Archive> void load(Archive& ar)
  {
    // TODO GET VIEW LABEL
    // Deserialize data into local variables and a vector
    ar(n_allocated_element, n_used_elements(), extra_allocation_factor);

    std::vector<T> data_vector;
    ar(data_vector); // Deserialize the data into the vector

    // if (n_used_elements > 0)
    {
      auto tmpdata = Kokkos::View<T*, Layout, HostSpace, Kokkos::MemoryTraits<Kokkos::Restrict>>(
          data_vector.data(), data_vector.size());
      _owned_data = Kokkos::create_mirror_view_and_copy(Space(), tmpdata);
      Kokkos::resize(_owned_data, n_allocated_element);
    }
  }

protected:
  void set_n_used_elements(uint64_t value)
  {
    KOKKOS_ASSERT(n_used_elements() < n_allocated_element);
    n_used_elements() = value;
  }

private:
  // public:
  uint64_t n_allocated_element = 0; ///< Number of allocated elements.
  // mutable uint64_t n_used_elements = 0; ///< Number of elements currently in use.

  Kokkos::View<uint64_t, Kokkos::SharedSpace> n_used_elements;

  double extra_allocation_factor =
      default_extra_allocation_factor; ///< Factor used for resizing the vector.

  static constexpr double default_extra_allocation_factor = 1.5;

  /**
   * @brief Allocates memory for the vector based on the new size.
   *
   * @param new_size The new size to allocate.
   * @return The size that was allocated.
   */
  size_t __allocate__(std::size_t new_size);
};

template <typename T, typename Space>
KokkosVector<T, Space>::KokkosVector(const size_t capacity, bool alloc, std::string label)
    : _owned_data(Kokkos::view_alloc(label + "_owned_data", Kokkos::WithoutInitializing)),
      n_used_elements("n_used_elements")
{
  __allocate__(capacity);
  if (alloc)
  {
    n_used_elements() = capacity;
  }
  else
  {
    n_used_elements() = 0;
  }
}

template <typename T, typename Space> void KokkosVector<T, Space>::clear() noexcept
{
  n_used_elements() = 0;
}

template <typename T, typename Space> void KokkosVector<T, Space>::reset()
{
  n_used_elements() = 0;
  n_allocated_element = 0;
  Kokkos::resize(_owned_data, n_allocated_element);
}

template <typename T, typename Space> void KokkosVector<T, Space>::resize(size_t n)
{
  this->n_used_elements() = 0;
  this->__allocate__(n);
}

template <typename T, typename Space>
KokkosVector<T, Space> KokkosVector<T, Space>::with_capacity(std::size_t capacity)
{
  auto rhs = KokkosVector(capacity);
  rhs.n_used_elements() = 0;
  return rhs;
}

// template <typename T, typename Space>
// KOKKOS_FUNCTION bool KokkosVector<T, Space>::emplace(T&& d) const
// {
//   const auto local_used_size = n_used_elements();
  
//   if (local_used_size < n_allocated_element)
//   {
//     Kokkos::atomic_increment(&n_used_elements());
//     _owned_data(local_used_size) = std::forward<T>(d);
//     return true;
//   }
//   return false;
// }


template <typename T, typename Space>
KOKKOS_FUNCTION bool KokkosVector<T, Space>::emplace(T&& d) const
{
  
  if (Kokkos::atomic_load(&n_used_elements()) < n_allocated_element)
  {
    _owned_data(Kokkos::atomic_fetch_inc(&n_used_elements())) = std::forward<T>(d);
    return true;
  }
  return false;
}




template <typename T, typename Space>
[[nodiscard]] double KokkosVector<T, Space>::get_allocation_factor() const noexcept
{
  return extra_allocation_factor;
}
template <typename T, typename Space>
void KokkosVector<T, Space>::set_allocation_factor(double value) noexcept
{
  KOKKOS_ASSERT(value >= 1);
  this->extra_allocation_factor = value;
}

template <typename T, typename Space>
void KokkosVector<T, Space>::insert(const KokkosVector<T, Space>& rhs)
{

  const auto original_size = n_used_elements();
  const auto n_add_item = rhs.size();
  __allocate__(original_size + rhs.size());

  auto data = this->data();
  
  Kokkos::parallel_for(
      "InsertNew",
      Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_add_item),
      KOKKOS_LAMBDA(const size_t i) { data(original_size + i) = rhs._owned_data(i); });
  Kokkos::fence();

  n_used_elements() += n_add_item;
}


template <typename T, typename Space>
size_t KokkosVector<T, Space>::__allocate__(const std::size_t new_size)
{
  if (new_size > 0)
  {
    if (new_size >= n_allocated_element)
    {
      const auto new_allocated_size = static_cast<std::size_t>(
          std::ceil(static_cast<double>(new_size) * extra_allocation_factor));

      n_allocated_element = new_allocated_size;
      Kokkos::resize(_owned_data, n_allocated_element);
    }
  }
  return new_size;
}

#endif