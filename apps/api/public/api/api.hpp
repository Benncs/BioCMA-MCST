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
  Handle(const Handle&) = delete;
  Handle(Handle&&) = default;
  Handle& operator=(const Handle&) = delete;
  Handle& operator=(Handle&&) = default;

  static std::optional<std::unique_ptr<Handle>>
  init(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_process);
  Handle() = default;
  ~Handle();

  bool regisiter_initial_condition();                          // TODO
  bool register_feed(Simulation::Feed::SimulationFeed&& feed); // TODO

  void register_parameters();
  int apply(bool to_load);
  void register_parameters(Core::UserControlParameters&& params);
  bool register_result_path(std::string_view path);
  bool register_cma_path(std::string_view path,bool recursive=false);
  bool register_serde(std::string_view path);

  bool register_model_name(std::string_view path);
  //flow, {glucose_c}, {0}, {0}, Simulation::Feed::Constant{}
  bool set_feed_constant(double _f, std::span<double> _target, std::span<std::size_t> _position, std::span<std::size_t> _species,bool gas=false);

  [[nodiscard]] int get_id() const;

  bool exec();

private:
  int id{};
  Handle(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_proces);
  Core::CaseData _data;
  Core::UserControlParameters params;
  bool loaded = false;
  bool applied = false;
  bool registered = false;
  std::optional<Simulation::Feed::SimulationFeed> feed=std::nullopt;
};

// } //namespace Api
#endif
