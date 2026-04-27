#ifndef __MC_EVENTS_HPP__
#define __MC_EVENTS_HPP__

#include "biocma_cst_config.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
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
    NewParticle = 0, ///< Spawn new particle
    Exit,            ///< Remove particle from list due to move in domain
    Move,            ///< Move in domain
    Death,           ///< Remove particle from list
    Overflow,
    ChangeWeight, ///< Update in weight
    __COUNT__
  };

  constexpr size_t number_event_type = static_cast<std::size_t>(
      EventType::__COUNT__); //< Number of different events
  /**
   * @brief inline getter, converts event to its value in order to be use as
   * array index
   */
  template <EventType event>
  KOKKOS_INLINE_FUNCTION consteval size_t
  event_index()
  {
    return static_cast<size_t>(event);
  }

  // class Events
  // {
  // public:
  //   using value_type = Kokkos::Array<std::size_t, number_event_type>;
  //   using view_type = Kokkos::View<value_type,
  //   Kokkos::DefaultExecutionSpace>;

  //   value_type counts;
  //   // std::size_t counts[number_event_type];

  //   KOKKOS_INLINE_FUNCTION void
  //   add_into(Events& dst) const
  //   {
  //     for (std::size_t i = 0; i < number_event_type; ++i)
  //     {
  //       dst.counts[i] += counts[i];
  //     }
  //   }

  //   KOKKOS_INLINE_FUNCTION
  //   Events() : counts()
  //   {
  //     init();
  //   }

  //   KOKKOS_INLINE_FUNCTION
  //   Events(const Events& rhs) : counts()
  //   {
  //     for (std::size_t i = 0; i < number_event_type; i++)
  //     {
  //       counts[i] = rhs.counts[i];
  //     }
  //   }

  //   template <EventType Event>
  //   KOKKOS_INLINE_FUNCTION void
  //   inc()
  //   {
  //     counts[event_index<Event>()]++;
  //   }
  //   template <EventType Event>
  //   KOKKOS_INLINE_FUNCTION void
  //   add(std::size_t val)
  //   {
  //     counts[event_index<Event>()] += val;
  //   }

  //   KOKKOS_INLINE_FUNCTION
  //   void
  //   init()
  //   {
  //     for (std::size_t i = 0; i < number_event_type; i++)
  //     {
  //       counts[i] = 0;
  //     }
  //   }
  // };

  // struct TallyReducer
  // {
  //   using reducer = TallyReducer;
  //   using value_type = Events;
  //   using result_view_type = Kokkos::View<value_type, Kokkos::HostSpace>;

  // private:
  //   result_view_type result_;

  // public:
  //   explicit TallyReducer() : result_("reducer")
  //   {
  //   }

  //   KOKKOS_INLINE_FUNCTION
  //   void
  //   join(value_type& dst, const value_type& src) const
  //   {
  //     src.add_into(dst);
  //   }

  //   KOKKOS_INLINE_FUNCTION
  //   void
  //   init(value_type& val) const
  //   {
  //     val.init();
  //   }

  //   [[nodiscard]] KOKKOS_INLINE_FUNCTION value_type&
  //   reference() const
  //   {
  //     return *result_.data();
  //   }

  //   [[nodiscard]] KOKKOS_INLINE_FUNCTION result_view_type
  //   view() const
  //   {
  //     return result_;
  //   }

  //   [[nodiscard]] KOKKOS_INLINE_FUNCTION bool
  //   references_scalar() const
  //   {
  //     return true;
  //   }
  // };

  /**
   * @brief Use to count events that occurs during Monte-Carlo processing cycles
   */
  struct EventContainer
  {

    // Use SharedHostPinnedSpace as this struct is meant to be small enough to
    // be shared between Host and Device According to SharedHostPinnedSpace
    // documentation, the size of this data can fit into one cache line so
    // transfer is not a botteneck
    Kokkos::View<std::size_t[number_event_type], Kokkos::SharedSpace> // NOLINT
        _events; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    /**
     * @brief Default container, initalise counter
     */
    EventContainer() : _events("events")
    {

      Kokkos::deep_copy(_events, 0); // Ensure all event to 0 occurence
    }
    /**
     * @brief Get std const view of _events counter
     */
    [[nodiscard]] std::span<std::size_t>
    get_span() const
    {
      // As we use SharedHostPinnedSpace, we can deal with classic std
      // containers to use host manipulation
      return { _events.data(), number_event_type };
    }

    /**
     * @brief Reset counters
     */
    KOKKOS_INLINE_FUNCTION void
    clear() const
    {
      Kokkos::deep_copy(_events, 0);
    }

    /**
     * @brief Getter to specific event counter
     * @tparam Event to get
     * @warning This function is thread_safe because we get a copy of the
     * current counter
     */
    template <EventType event>
    [[nodiscard]] constexpr std::size_t
    get() const
    {
      static_assert(event != EventType::__COUNT__,
                    "Count is not a valid event");
      return _events(event_index<event>());
    }

    /**
     * @brief Increment specific event counter
     *
     * @warning This function is thread_safe and can be use in kernel
     * @tparam Event to get
     *
     */
    template <EventType event>
    KOKKOS_FORCEINLINE_FUNCTION constexpr void
    incr() const
    {
      static_assert(event != EventType::__COUNT__,
                    "Count is not a valid event");
      Kokkos::atomic_add(&_events[event_index<event>()], 1);
    }

    template <EventType event>
    KOKKOS_FORCEINLINE_FUNCTION constexpr void
    add(size_t val) const
    {
      static_assert(event != EventType::__COUNT__,
                    "Count is not a valid event");
      Kokkos::atomic_add(&_events[event_index<event>()], val);
    }

    template <EventType event>
    KOKKOS_FORCEINLINE_FUNCTION constexpr void
    wrap_incr() const
    {
      if constexpr (AutoGenerated::FlagCompileTime::enable_event_counter)
      {
        incr<event>();
      }
    }

    template <class Archive>
    void
    save(Archive& ar) const
    {
      std::array<std::size_t, number_event_type> array{};

      auto rd = std::span<std::size_t>(_events.data(), number_event_type);

      std::copy(rd.begin(), rd.end(), array.begin());
      assert(rd[0] == _events[0] && rd[0] == array[0]);

      ar(array);
    }

    template <class Archive>
    void
    load(Archive& ar)
    {

      std::array<std::size_t, number_event_type> array{};
      ar(array);

      auto rd = std::span<std::size_t>(_events.data(), number_event_type);

      std::copy(array.begin(), array.end(), rd.begin());
      assert(rd[0] == _events[0]);
    }
  };

}; // namespace MC

#endif //__MC_EVENTS_HPP__
