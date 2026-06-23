#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>
#include <simulation/feed_descriptor.hpp>

constexpr double TEST_TOLERANCE = 1e-12;

void
basic_insertion()
{
  auto constant_feed = Simulation::Feed::FeedDescriptor{
    1, { { 1, 0 } }, 1, 1, Simulation::Feed::Constant{}
  };

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

void
test_update_constant()
{
  double original_flow = 1;
  auto constant_feed = Simulation::Feed::FeedDescriptor{
    original_flow, { { 1, 0 } }, 1, 1, Simulation::Feed::Constant{}
  };

  double t = 0;

  for (int i = 0; i < 10; ++i)
  {
    constant_feed.update(t, 0.1);
    assert(constant_feed.flow == original_flow);

    t += 0.1;
  }
}

void
test_update_linear()
{
  double original_flow = 1;
  double df = 0.1;
  auto lin = Simulation::Feed::Linear{ .f0 = 1, .df = df };
  double n_it = 11;
  double d_t = 0.1;
  double t_end = n_it * d_t;

  auto feed = Simulation::Feed::FeedDescriptor{
    original_flow, { { 1, 0 } }, 1, 1, lin
  };

  double t = 0;

  for (int i = 0; i < n_it + 1; ++i)
  {
    auto expected = original_flow + t * df;

    feed.update(t, d_t);

    assert(std::abs(feed.flow - expected) < TEST_TOLERANCE);
    t += d_t;
  }
  auto expected = original_flow + t_end * df;
  std::cerr << feed.flow << std::endl;
  std::cerr << expected << std::endl;
  assert(std::abs(feed.flow - expected) < TEST_TOLERANCE);
}

void
test_update_exp()
{
  double original_flow = 1;
  double alpha = 0.1;
  auto exp = Simulation::Feed::Exponential{ .f0 = 1, .alpha = alpha };
  double n_it = 11;
  double d_t = 0.1;
  double t_end = n_it * d_t;

  auto feed = Simulation::Feed::FeedDescriptor{
    original_flow, { { 1, 0 } }, 1, 1, exp
  };

  double t = 0;

  for (int i = 0; i < n_it + 1; ++i)
  {
    auto expected = original_flow + std::exp(t * alpha);

    feed.update(t, d_t);

    assert(std::abs(feed.flow - expected) < TEST_TOLERANCE);
    t += d_t;
  }
  auto expected = original_flow + std::exp(t_end * alpha);
  std::cerr << feed.flow << std::endl;
  std::cerr << expected << std::endl;
  assert(std::abs(feed.flow - expected) < TEST_TOLERANCE);
}

int
main()
{

  basic_insertion();
  test_update_constant();
  test_update_linear();
  test_update_exp();
  return 0;
}