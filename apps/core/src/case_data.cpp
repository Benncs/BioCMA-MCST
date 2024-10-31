#include "biocma_cst_config.hpp"
#include "common/execinfo.hpp"

#include "core/global_initaliser.hpp"
#include "post_process.hpp"
#include <core/case_data.hpp>
#include <host_specific.hpp>
#include <optional>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <worker_specific.hpp>

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
    filename << params.user_params.results_file_name << "_partial_" << exec.current_rank << ".h5";
    return filename.str();
  }

  void init_partial_exporter(CORE_DE::PartialExporter &partial_exporter, const CaseData &case_data)
  {

    const auto [_, n_compartment] = case_data.simulation->getDimensions();

    partial_exporter.init_fields(case_data.params.user_params.number_exported_result, n_compartment);

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

    std::cout << " prima  " << case_data.simulation->mc_unit->domain[0].n_cells << std::endl;
    std::cout << " prima  " << case_data.simulation->mc_unit->init_weight << std::endl;
    CORE_DE::PartialExporter partial_exporter(ExecInfo(), filename);
    std::cout << " poi  " << case_data.simulation->mc_unit->domain[0].n_cells << std::endl;
    std::cout << " poi  " << case_data.simulation->mc_unit->init_weight << std::endl;
    init_partial_exporter(partial_exporter, case_data);

    auto *const sim = case_data.simulation.get();
    f_run(case_data.exec_info, std::move(*sim), case_data.params, std::move(case_data.transitioner), partial_exporter);

    do_serde(case_data);
  }

  std::optional<Core::CaseData> load(const ExecInfo &exec, SimulationParameters &params)
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

    bool ok_init = SerDe::load_simulation(
        gi, case_data, "/mnt/c/Users/casale/Documents/code/cpp/kokkos_biomc/results/debug_serde/debug_serde_serde_0");

    if (!gi.check_init_terminate() || !ok_init)
    {
      return std::nullopt;
    }

    case_data.params = params;
    // TODO: Check integreity between transitionner and read value (transitionner and simulation tpf + n_species)
    std::cout << " load:   " << case_data.simulation->mc_unit->domain[0].n_cells << std::endl;
    return case_data;
  }

} // namespace Core