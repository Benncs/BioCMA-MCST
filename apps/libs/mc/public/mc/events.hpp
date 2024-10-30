#ifndef __MC_EVENTS_HPP__
#define __MC_EVENTS_HPP__

#include "common/kokkos_vector.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Macros.hpp>
#include <common/execinfo.hpp>
#include <cstddef>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <span>

namespace MC
{
  /**
   * @brief Enumeration that represents events that can occurs during a
   * Monte-Carlo cycle
   */
  enum class EventType : char
  {
    ChangeWeight = 0, ///< Update in weight
    NewParticle = 1,  ///< Spawn new particle
    Death = 2,        ///< Remove particle from list
    Move = 3,         ///< Move in domain
    Exit = 4          ///< Remove particle from list due to move in domain
  };

  constexpr size_t number_event_type = 5; //< Number of different events

  /**
   * @brief inline getter, converts event to its value in order to be use as
   * array index
   */
  template <EventType event> KOKKOS_INLINE_FUNCTION consteval size_t event_index()
  {
    return static_cast<size_t>(event);
  }

  /**
   * @brief Use to count events that occurs during Monte-Carlo processing cycles
   */
  struct alignas(ExecInfo::cache_line_size) EventContainer
  {

    // Use SharedHostPinnedSpace as this struct is meant to be small enough to
    // be shared between Host and Device According to SharedHostPinnedSpace
    // documentation, the size of this data can fit into one cache line so
    // transfer is not a botteneck
    Kokkos::View<size_t[number_event_type], Kokkos::SharedHostPinnedSpace>
        _events; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    static_assert(sizeof(size_t[number_event_type]) <= ExecInfo::cache_line_size, "Size of Eventcontainer::_event fits into cache line");

    /**
     * @brief Default container, initalise counter
     */
    EventContainer()
    {
      _events = Kokkos::View<size_t[number_event_type],
                             Kokkos::SharedHostPinnedSpace>("events"); // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
      Kokkos::deep_copy(_events, 0);                                   // Ensure all event to 0 occurence
    }
    /**
     * @brief Get std const view of _events counter
     */
    [[nodiscard]] std::span<size_t> get_span() const
    {
      // As we use SharedHostPinnedSpace, we can deal with classic std
      // containers to use host manipulation
      return {_events.data(), number_event_type};
    }

    /**
     * @brief Reset counters
     */
    KOKKOS_INLINE_FUNCTION void clear() const
    {
      Kokkos::deep_copy(_events, 0);
    }

    /**
     * @brief Getter to specific event counter
     * @tparam Event to get
     * @warning This function is thread_safe because we get a copy of the
     * current counter
     */
    template <EventType event> [[nodiscard]] constexpr size_t get() const
    {
      return _events(event_index<event>());
    }

    /**
     * @brief Getter to specific mutable event counter
     * @tparam Event to get
     * @warning This function is NOT thread_safe because we get a raw reference
     * on the current counter
     */
    template <EventType event> [[nodiscard]] constexpr size_t &get_mut()
    {
      return _events(event_index<event>());
    }

    /**
     * @brief Increment specific event counter
     *
     * @warning This function is thread_safe and can be use in kernel
     * @tparam Event to get
     *
     */
    template <EventType event> KOKKOS_INLINE_FUNCTION constexpr void incr() const
    {
      Kokkos::atomic_increment(&_events(event_index<event>()));
    }

    /**
     * @brief Transform a linear contiguous counter data obtained via MPI gather
     * into EventContainer object
     * @param _data obtained via multiple EventContainer merged together
     * @warning _data size has to be a multiple of number_event_type
     */
    static EventContainer reduce(std::span<size_t> _data);

    // Either for save and load, as number_event_type is small and data is located in host, we can loop over array

    template <class Archive> void save(Archive &ar) const
    {
      std::array<size_t, number_event_type> array{};
      for (size_t i = 0; i < number_event_type; ++i)
      {
        array[i] = _events[i];
      }

      ar(array);
    }

    template <class Archive> void load(Archive &ar)
    {
      std::array<size_t, number_event_type> array{};
      ar(array);
      for (size_t i = 0; i < number_event_type; ++i)
      {
        _events(i) = array[i];
      }
    }
  };

  static_assert(sizeof(EventContainer) <= 2 * ExecInfo::cache_line_size, "Check size of Eventcontainer");

}; // namespace MC

#endif //__MC_EVENTS_HPP__