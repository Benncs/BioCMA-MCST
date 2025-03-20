

#ifndef __CORE__SERDE_HPP__
#define __CORE__SERDE_HPP__
#include <core/global_initaliser.hpp>

#include <core/case_data.hpp>
#include <string_view>
#ifdef USE_CEAREAL
namespace SerDe
{
  void save_simulation(const Core::CaseData& case_data);
  bool load_simulation(Core::GlobalInitialiser& gi,
                       Core::CaseData& case_data,
                       std::string_view ser_filename);
} // namespace SerDe

#endif // USE_CEAREAL

#ifdef USE_CEAREAL
#  include <cereal/archives/binary.hpp>
// template <typename T = void>
// std::enable_if_t<AutoGenerated::FlagCompileTime::use_cereal_serde, T>
// do_serde(Core::CaseData& case_data)
// {
//   std::cout << "Serialization..." << std::endl;
//   SerDe::save_simulation(case_data);
// }
template <typename T = void>
  requires AutoGenerated::FlagCompileTime::use_cereal_serde
void do_serde(Core::CaseData& case_data)
{
  std::cout << "Serialization..." << std::endl;
  SerDe::save_simulation(case_data);
}

template <typename T = void>
  requires AutoGenerated::FlagCompileTime::use_cereal_serde
std::optional<Core::CaseData> impl_load(const ExecInfo& exec,
                                        const Core::UserControlParameters&& params,
                                        std::optional<Simulation::Feed::SimulationFeed> feed)
{
  Core::CaseData case_data;
  case_data.exec_info = exec;
  // if (exec.n_rank > 1)
  // {
  //   throw std::runtime_error("Serde for MPI Not implemented yet");
  // }

  Core::GlobalInitialiser gi(exec, params);
  auto transition = gi.init_transitionner();
  if (!transition.has_value())
  {
    return std::nullopt;
  }

  case_data.transitioner = std::move(*transition);
  gi.init_feed(std::move(feed));
  try
  {
    if (!params.serde_file.has_value())
    {
      return std::nullopt;
    }
    std::string serde_filename = *params.serde_file + std::to_string(exec.current_rank) + ".raw";
    std::cout << serde_filename << std::endl;
    const bool ok_init = SerDe::load_simulation(gi, case_data, serde_filename);

    if (!gi.check_init_terminate() || !ok_init)
    {
      return std::nullopt;
    }
  }
  catch (std::exception& e)
  {
    auto err = "CORE::load::load_simulation:" + std::string(e.what());
    throw std::runtime_error(err);
  }

  case_data.params = gi.get_parameters();
  // TODO: Check integreity between transitionner and read value (transitionner and simulation tpf
  // + n_species)

  return case_data;
}

#endif

template <typename T = void>
  requires(!AutoGenerated::FlagCompileTime::use_cereal_serde)
void do_serde(Core::CaseData& /*arg*/)
{
  // NOP
}

template <typename T = void>
  requires(!AutoGenerated::FlagCompileTime::use_cereal_serde)
[[maybe_unused]] std::optional<Core::CaseData>
impl_load([[maybe_unused]] const ExecInfo& exec,
          [[maybe_unused]] const Core::UserControlParameters&& params,
          [[maybe_unused]] std::optional<Simulation::Feed::SimulationFeed> feed)
{
  return std::nullopt;
}

#endif
