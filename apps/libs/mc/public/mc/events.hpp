#ifndef __MC_EVENTS_HPP__
#define __MC_EVENTS_HPP__

#include "Kokkos_Atomic.hpp"
#include "Kokkos_Core_fwd.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <common/execinfo.hpp>
#include <cstddef>
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
    Kokkos::View<std::size_t[number_event_type], Kokkos::SharedSpace> // NOLINT
        _events; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    // NOLINTBEGIN
    static_assert(sizeof(std::size_t[number_event_type]) <= ExecInfo::cache_line_size,
                  "Size of Eventcontainer::_event fits into cache line");
    // NOLINTEND
    /**
     * @brief Default container, initalise counter
     */
    EventContainer():_events("events")
    {
  
      Kokkos::deep_copy(_events, 0); // Ensure all event to 0 occurence
    }
    /**
     * @brief Get std const view of _events counter
     */
    [[nodiscard]] std::span<std::size_t> get_span() const
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
    template <EventType event> [[nodiscard]] constexpr std::size_t get() const
    {
      return _events(event_index<event>());
    }

    /**
     * @brief Getter to specific mutable event counter
     * @tparam Event to get
     * @warning This function is NOT thread_safe because we get a raw reference
     * on the current counter
     */
    template <EventType event> [[nodiscard]] constexpr std::size_t& get_mut()
    {
      return _events[event_index<event>()];
    }

    /**
     * @brief Increment specific event counter
     *
     * @warning This function is thread_safe and can be use in kernel
     * @tparam Event to get
     *
     */
    template <EventType event> KOKKOS_FORCEINLINE_FUNCTION constexpr void incr() const
    {
      Kokkos::atomic_add(&_events[event_index<event>()], 1);
    }

    template <EventType event> KOKKOS_FORCEINLINE_FUNCTION constexpr void add(size_t val) const
    {
      Kokkos::atomic_add(&_events[event_index<event>()], val);
    }

    /**
     * @brief Transform a linear contiguous counter data obtained via MPI gather
     * into EventContainer object
     * @param _data obtained via multiple EventContainer merged together
     * @warning _data size has to be a multiple of number_event_type
     */
    static EventContainer reduce(std::span<std::size_t> _data);

    template <class Archive> void save(Archive& ar) const
    {
      std::array<std::size_t, number_event_type> array{};

      auto rd = std::span<std::size_t>(_events.data(), number_event_type);

      std::copy(rd.begin(), rd.end(), array.begin());
      assert(rd[0] == _events[0] && rd[0] == array[0]);

      ar(array);
    }

    template <class Archive> void load(Archive& ar)
    {

      std::array<std::size_t, number_event_type> array{};
      ar(array);

      auto rd = std::span<std::size_t>(_events.data(), number_event_type);

      std::copy(array.begin(), array.end(), rd.begin());
      assert(rd[0] == _events[0]);
    }
  };

  static_assert(sizeof(EventContainer) <= 2 * ExecInfo::cache_line_size,
                "Check size of Eventcontainer");

}; // namespace MC

#endif //__MC_EVENTS_HPP__