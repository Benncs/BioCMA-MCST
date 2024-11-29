#ifndef __BIOMC_API_HPP__
#define __BIOMC_API_HPP__

#include <core/case_data.hpp>
#include <core/simulation_parameters.hpp>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <span>
#include <string_view>

#include <api/results.hpp>
namespace Api
{
class SimulationInstance
{
public:
  SimulationInstance(const SimulationInstance&) = delete;
  SimulationInstance(SimulationInstance&&) = default;
  SimulationInstance& operator=(const SimulationInstance&) = delete;
  SimulationInstance& operator=(SimulationInstance&&) = default;

  // TODO Enable if def USE_MPI
  static std::optional<std::unique_ptr<SimulationInstance>>
  init(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_process) noexcept;

  static std::optional<std::unique_ptr<SimulationInstance>> init(uint64_t id,
                                                     uint32_t thread_per_process) noexcept;
  SimulationInstance() = default;
  ~SimulationInstance() = default;

  ApiResult apply(bool to_load) noexcept;
  ApiResult register_parameters(Core::UserControlParameters&& params) noexcept;
  bool register_result_path(std::string_view path);
  bool register_cma_path(std::string_view path, bool recursive = false);
  bool register_serde(std::string_view path);

  ApiResult register_model_name(std::string_view path);
  bool set_feed_constant(double _f,
                         std::span<double> _target,
                         std::span<std::size_t> _position,
                         std::span<std::size_t> _species,
                         bool gas = false);

  [[nodiscard]] int get_id() const;

  ApiResult exec() noexcept;

private:
  int id{};
  SimulationInstance(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_proces);
  Core::CaseData _data;
  Core::UserControlParameters params;
  bool loaded = false;
  bool applied = false;
  bool registered = false;
  std::optional<Simulation::Feed::SimulationFeed> feed = std::nullopt;
};

} //namespace Api
#endif
