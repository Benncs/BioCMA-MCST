#include <cassert>
#include <optional>
#include <simulation/feed_descriptor.hpp>

void basic_insertion()
{
  auto constant_feed =
      Simulation::Feed::FeedDescriptor{1, 1, 0, 1, 1, Simulation::Feed::Constant{}};

  auto feeds = Simulation::Feed::SimulationFeed();

  assert(feeds.liquid_feeds().size() == 0);

  feeds.add_liquid(std::move(constant_feed));

  {
    assert(feeds.liquid_feeds().size() == 1);
  }

  feeds.add_feed(std::move(constant_feed), Phase::Liquid);
  {
    assert(feeds.liquid_feeds().size() == 2);
  }

  feeds.add_feed(std::move(constant_feed), Phase::Gas);
  {
    assert(feeds.liquid_feeds().size() == 2);
    assert(feeds.gas_feeds().size() == 1);
  }

  feeds.add_gas(std::move(constant_feed));
  {
    assert(feeds.liquid_feeds().size() == 2);
    assert(feeds.gas_feeds().size() == 2);
  }
}

int main()
{

  basic_insertion();
  return 0;
}