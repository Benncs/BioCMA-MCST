#include "common/execinfo.hpp"


#include "post_process.hpp"
#include <core/case_data.hpp>
#include <host_specific.hpp>
#include <sstream>
#include <string>
#include <string_view>
#include <worker_specific.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>

namespace Core
{

  std::string exporter_filename(const ExecInfo &exec,
                                const SimulationParameters &params)
  {
    std::stringstream filename;
    filename 
             << params.user_params.results_file_name << "_partial_"
             << exec.current_rank << ".h5";
    return filename.str();
  }

  void init_partial_exporter(CORE_DE::PartialExporter &partial_exporter,
                             const CaseData &case_data)
  {

    const auto [_, n_compartment] = case_data.simulation->getDim();

    partial_exporter.init_fields(
        case_data.params.user_params.number_exported_result, n_compartment);

    auto probes = Simulation::Probes(0, 1, {"exit_time"});
    case_data.simulation->set_probes(std::move(probes));

    partial_exporter.write_number_particle(
        case_data.simulation->mc_unit->domain.getRepartition());

    PostProcessing::save_initial_particle_state(*case_data.simulation,
                                                partial_exporter);
  }

  void exec(CaseData &&case_data)
  {
#ifdef NO_MPI
    const auto f_run = host_process;
#else
    const auto f_run = (case_data.exec_info.current_rank == 0)
                           ? &host_process
                           : &workers_process;
#endif

    const auto filename =
        exporter_filename(case_data.exec_info, case_data.params);
    CORE_DE::PartialExporter partial_exporter(case_data.exec_info, filename);

    init_partial_exporter(partial_exporter, case_data);
    auto *const sim = case_data.simulation.release();
    f_run(case_data.exec_info,
          std::move(*sim),
          case_data.params,
          std::move(case_data.transitioner),
          partial_exporter);
  }

} // namespace Core