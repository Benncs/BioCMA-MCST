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
#include <string>
#include <string_view>
#include <variant>
// namespace Api
//
// {
struct Success
{
};
template <typename T> struct Result : protected std::variant<Success, T>
{
  explicit constexpr Result() noexcept : std::variant<Success, T>{Success{}} {};
  constexpr explicit Result(T const&& t) noexcept : std::variant<Success, T>{t}
  {
  }
  constexpr explicit operator bool() const noexcept
  {
    return valid();
  }

  [[nodiscard]] constexpr bool valid() const noexcept
  {
    return std::holds_alternative<Success>(*this);
  }
  [[nodiscard]] constexpr bool invalid() const noexcept
  {
    return !valid();
  }
  [[nodiscard]] constexpr auto get() const noexcept -> T
  {
    return (invalid() ? std::get<T>(*this) : T());
  }
};

struct ApiResult : Result<std::string>
{
  explicit ApiResult(std::string_view t) noexcept : Result<std::string>(std::string(t))
  {
  }

  explicit constexpr ApiResult() noexcept = default;
  constexpr int to_c_ret_code()
  {
    return (valid()) ? 0 : -1;
  }
};

class Handle
{
public:
  Handle(const Handle&) = delete;
  Handle(Handle&&) = default;
  Handle& operator=(const Handle&) = delete;
  Handle& operator=(Handle&&) = default;

  // TODO Enable if def USE_MPI
  static std::optional<std::unique_ptr<Handle>>
  init(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_process) noexcept;
  static std::optional<std::unique_ptr<Handle>> init(uint64_t id,
                                                     uint32_t thread_per_process) noexcept;
  Handle() = default;
  ~Handle() = default;

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
  Handle(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_proces);
  Core::CaseData _data;
  Core::UserControlParameters params;
  bool loaded = false;
  bool applied = false;
  bool registered = false;
  std::optional<Simulation::Feed::SimulationFeed> feed = std::nullopt;
};

// } //namespace Api
#endif
