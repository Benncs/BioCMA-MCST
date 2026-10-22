#include <Kokkos_Core.hpp>
#include <array>
#include <cassert>
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cstddef>
#include <cstdio>
#include <mc/events.hpp>
#include <sstream>

namespace
{
  constexpr std::size_t n_iter = 100000;

  // One distinct increment count per event so a cross-talking counter shows up
  // as a wrong total rather than a coincidentally matching one.
  constexpr std::array<std::size_t, MC::number_event_type> multiplier
      = { 1, 2, 3, 4, 5, 6 };

  template <MC::EventType Event, std::size_t Mult>
  void
  hammer(const MC::EventContainer& events)
  {
    Kokkos::parallel_for(
        "hammer",
        Kokkos::RangePolicy<>(0, n_iter),
        KOKKOS_LAMBDA(const std::size_t) {
          for (std::size_t k = 0; k < Mult; ++k)
          {
            events.incr<Event>();
          }
        });
    Kokkos::fence();
  }

  int failures = 0;

  void
  check(const char* what, std::size_t got, std::size_t expected)
  {
    if (got != expected)
    {
      std::printf("FAIL %-34s got %zu expected %zu\n", what, got, expected);
      ++failures;
    }
    else
    {
      std::printf("ok   %-34s %zu\n", what, got);
    }
  }
} // namespace

int
main()
{
  Kokkos::initialize();
  {
    MC::EventContainer events;

    hammer<MC::EventType::NewParticle, multiplier[0]>(events);
    hammer<MC::EventType::Exit, multiplier[1]>(events);
    hammer<MC::EventType::Move, multiplier[2]>(events);
    hammer<MC::EventType::Death, multiplier[3]>(events);
    hammer<MC::EventType::Overflow, multiplier[4]>(events);
    hammer<MC::EventType::ChangeWeight, multiplier[5]>(events);

    check("get<NewParticle>",
          events.get<MC::EventType::NewParticle>(),
          n_iter * multiplier[0]);
    check(
        "get<Exit>", events.get<MC::EventType::Exit>(), n_iter * multiplier[1]);
    check(
        "get<Move>", events.get<MC::EventType::Move>(), n_iter * multiplier[2]);
    check("get<Death>",
          events.get<MC::EventType::Death>(),
          n_iter * multiplier[3]);
    check("get<Overflow>",
          events.get<MC::EventType::Overflow>(),
          n_iter * multiplier[4]);
    check("get<ChangeWeight>",
          events.get<MC::EventType::ChangeWeight>(),
          n_iter * multiplier[5]);

    // get_span() must expose the same values, packed and in enum order.
    const auto packed = events.get_span();
    for (std::size_t i = 0; i < MC::number_event_type; ++i)
    {
      check("get_span()", packed[i], n_iter * multiplier[i]);
    }

    // add() must agree with incr() on the same padded slot.
    events.add<MC::EventType::Exit>(7);
    check("add<Exit>",
          events.get<MC::EventType::Exit>(),
          n_iter * multiplier[1] + 7);
    check("add<Exit> leaves Move",
          events.get<MC::EventType::Move>(),
          n_iter * multiplier[2]);

    const auto* base = events._events.data();
    const auto delta = static_cast<std::size_t>(&events._events(1) - base)
                       * sizeof(std::size_t);
    check("bytes between counters", delta, 64);

    // Round-trip: on-disk payload stays a packed array of number_event_type.
    std::ostringstream oss;
    {
      cereal::BinaryOutputArchive ar(oss);
      ar(events);
    }
    MC::EventContainer reloaded;
    {
      std::istringstream iss(oss.str());
      cereal::BinaryInputArchive ar(iss);
      ar(reloaded);
    }
    for (std::size_t i = 0; i < MC::number_event_type; ++i)
    {
      check("serde round-trip",
            reloaded.get_span()[i],
            packed[i] + (i == 1 ? 7 : 0));
    }

    check("payload bytes",
          oss.str().size(),
          sizeof(std::size_t) * MC::number_event_type);

    events.clear();
    for (const auto value : events.get_span())
    {
      check("clear()", value, 0);
    }
  }
  Kokkos::finalize();

  std::printf(failures == 0 ? "\nALL PASS\n" : "\n%d FAILURE(S)\n", failures);
  return failures == 0 ? 0 : 1;
}
