#ifndef __BIOMC_API_HPP__
#define __BIOMC_API_HPP__

#include "core/simulation_parameters.hpp"
#include "simulation/feed_descriptor.hpp"
#include <core/case_data.hpp>
#include <cstdlib>
#include <memory>
#include <optional>
// namespace Api
// {

class Handle
{
public:
  static std::optional<std::unique_ptr<Handle>> init(uint32_t n_rank, uint32_t current_rank, uint64_t id,uint32_t thread_per_process);
  static std::optional<std::unique_ptr<Handle>> load_mock(uint32_t n_rank, uint32_t current_rank);
  Handle() = default;
  ~Handle();
  void register_parameters();
  void apply(bool to_load);

  void register_parameters(Core::UserControlParameters &&params);
  bool regisiter_initial_condition();                          // TODO
  bool register_feed(Simulation::Feed::SimulationFeed &&feed); // TODO
  bool register_result_path(std::string_view path);
  bool register_cma_path(std::string_view path);
  [[nodiscard]] int get_id() const
  {
    return id;
  }

  bool exec();

private:
  int id{};
  Handle(uint32_t n_rank, uint32_t current_rank, uint64_t id,uint32_t thread_per_proces);
  Core::CaseData _data;
  Core::UserControlParameters params;
  bool loaded = false;
  bool applied = false;
  bool registered = false;
};

// } //namespace Api
#endif
