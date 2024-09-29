#include "post_process.hpp"
#include <core/case_data.hpp>
#include <host_specific.hpp>
#include <string>
#include <string_view>
#include <worker_specific.hpp>

namespace Core
{

  void exec(CaseData &&case_data)
  {
#ifdef NO_MPI
    const auto f_run = host_process;
#else
    const auto f_run = (case_data.exec_info.current_rank == 0)
                           ? &host_process
                           : &workers_process;
#endif

    //TODO Put this elsewere, clean function for either main and partial name
    auto filename =
        "./results/" + case_data.params.user_params.results_file_name +
        "_partial_" + std::to_string(case_data.exec_info.current_rank) + ".h5";
    CORE_DE::PartialExporter partial_exporter(case_data.exec_info, filename);
    const auto [n_compartment, n_species] = case_data.simulation->getDim();

    partial_exporter.init_fields(
        case_data.params.user_params.number_exported_result, n_compartment);

    auto probes = Simulation::Probes(0, 1, {"exit_time"});
    case_data.simulation->set_probes(std::move(probes));

    auto *const sim = case_data.simulation.release();

    partial_exporter.write_number_particle(
        sim->mc_unit->domain.getRepartition());

    PostProcessing::save_initial_particle_state(*sim, partial_exporter);

    f_run(case_data.exec_info,
          std::move(*sim),
          case_data.params,
          std::move(case_data.transitioner),
          partial_exporter);
  }

} // namespace Core