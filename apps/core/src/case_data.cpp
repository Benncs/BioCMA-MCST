#include "biocma_cst_config.hpp"
#include "common/execinfo.hpp"

#include "core/global_initaliser.hpp"
#include "post_process.hpp"
#include "simulation/feed_descriptor.hpp"
#include <core/case_data.hpp>
#include <host_specific.hpp>
#include <optional>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <worker_specific.hpp>

#ifdef NO_MPI
namespace MPI_W
{
  bool is_initialized() noexcept
  {
    return false;
  }
} // namespace MPI_W
#else
#  include "mpi_w/message_t.hpp"
#endif

#ifdef USE_CEAREAL
#  include "serde.hpp"
#  include <cereal/archives/binary.hpp>
template <typename T = void>
std::enable_if_t<AutoGenerated::FlagCompileTime::use_cereal_serde, T> do_serde(Core::CaseData &case_data)
{
  SerDe::save_simulation(case_data);
}

#endif

template <typename T = void>
std::enable_if_t<!AutoGenerated::FlagCompileTime::use_cereal_serde, T> do_serde(Core::CaseData & /**/)
{
  // Nop
}

namespace Core
{

  std::string exporter_filename(const ExecInfo &exec, const SimulationParameters &params)
  {
    std::stringstream filename;
    filename << params.results_file_name << "_partial_" << exec.current_rank << ".h5";
    return filename.str();
  }

  void init_partial_exporter(CORE_DE::PartialExporter &partial_exporter, const CaseData &case_data)
  {

    const auto [_, n_compartment] = case_data.simulation->getDimensions();

    partial_exporter.init_fields(case_data.params.number_exported_result, n_compartment);

    auto probes = Simulation::Probes(0, 1, {"exit_time"});
    case_data.simulation->set_probes(std::move(probes));

    partial_exporter.write_number_particle(case_data.simulation->mc_unit->domain.getRepartition());

    PostProcessing::save_initial_particle_state(*case_data.simulation, partial_exporter);
  }

  void exec(CaseData &&case_data)
  {
#ifdef NO_MPI
    const auto f_run = host_process;
#else
    const auto f_run = (case_data.exec_info.current_rank == 0) ? &host_process : &workers_process;
#endif

    const auto filename = exporter_filename(case_data.exec_info, case_data.params);

    CORE_DE::PartialExporter partial_exporter(ExecInfo(), filename);

    init_partial_exporter(partial_exporter, case_data);

    auto *const sim = case_data.simulation.get();
    f_run(case_data.exec_info, std::move(*sim), case_data.params, std::move(case_data.transitioner), partial_exporter);

    do_serde(case_data);
  }

  std::optional<Core::CaseData> load(const ExecInfo &exec, const UserControlParameters &&params,std::optional<Simulation::Feed::SimulationFeed> feed)
  {
    CaseData case_data;
    case_data.exec_info = exec;
    if (exec.n_rank > 1)
    {
      throw std::runtime_error("Serde for MPI Not implemented yet");
    }

    GlobalInitialiser gi(exec, params);
    auto transition = gi.init_transitionner();
    if (!transition.has_value())
    {
      return std::nullopt;
    }

    case_data.transitioner = std::move(*transition);
    gi.init_feed(std::move(feed));
    try
    {
      const bool ok_init = SerDe::load_simulation(gi, case_data, *params.serde_file);

      if (!gi.check_init_terminate() || !ok_init)
      {
        return std::nullopt;
      }
    }
    catch (std::exception &e)
    {
      auto err = "CORE::load::load_simulation:" + std::string(e.what());
      throw std::runtime_error(err);
    }

    case_data.params = gi.get_parameters();
    // TODO: Check integreity between transitionner and read value (transitionner and simulation tpf + n_species)

    return case_data;
  }

} // namespace Core
