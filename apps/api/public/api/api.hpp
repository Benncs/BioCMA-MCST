#ifndef __BIOMC_API_HPP__
#define __BIOMC_API_HPP__

#include <core/case_data.hpp>
#include "core/simulation_parameters.hpp"
#include "simulation/feed_descriptor.hpp"
#include <cstdlib>
#include <memory>
#include <optional>
// namespace Api
// {
class Handle
{
public:
  static std::optional<std::unique_ptr<Handle>> init(uint32_t n_rank,uint32_t current_rank);
  static std::optional<std::unique_ptr<Handle>> load(uint32_t n_rank,uint32_t current_rank);
  Handle() = default;
  ~Handle();
  void register_parameters();
  void apply();
  // void register_parameters(Core::SimulationParameters &&params);
  bool regisiter_initial_condition();
  bool register_feed(Simulation::Feed::SimulationFeed &&feed);
  bool exec();
  int id;

private:
  Core::CaseData _data;
  bool loaded = false;
  bool applied = false;
  bool registered=false;
};
// } //namespace Api
#endif
