#include <biocma_cst_config.hpp>
#include <common/execinfo.hpp>
#include <core/case_data.hpp>
#include <core/global_initaliser.hpp>
#include <host_specific.hpp>
#include <optional>
#include <impl_post_process.hpp>
#include <serde.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <worker_specific.hpp>

namespace Core
{

  std::string exporter_filename(const ExecInfo& exec, const SimulationParameters& params)
  {
    std::stringstream filename;
    filename << params.results_file_name << "_partial_" << exec.current_rank << ".h5";
    return filename.str();
  }

  void init_partial_exporter(Core::PartialExporter& partial_exporter, const CaseData& case_data)
  {
    const auto [_, n_compartment] = case_data.simulation->getDimensions();
    partial_exporter.init_fields(case_data.params.number_exported_result, n_compartment);
    auto probes = Simulation::ProbeAutogeneratedBuffer();
    case_data.simulation->set_probes(std::move(probes));
    partial_exporter.write_number_particle(case_data.simulation->mc_unit->getRepartition());
    PostProcessing::save_particle_state(*case_data.simulation, partial_exporter);

  }

  void exec(CaseData&& case_data)
  {
#ifdef NO_MPI
    const auto f_run = host_process;
#else
    const auto f_run = (case_data.exec_info.current_rank == 0) ? &host_process : &workers_process;
#endif

    const auto filename = exporter_filename(case_data.exec_info, case_data.params);
    Core::PartialExporter partial_exporter(ExecInfo(), filename);
    init_partial_exporter(partial_exporter, case_data);
    auto* const sim = case_data.simulation.get();
     
    f_run(case_data.exec_info,
          *sim,
          case_data.params,
          std::move(case_data.transitioner),
          partial_exporter);


    if constexpr (AutoGenerated::FlagCompileTime::use_probe)
    {
      PostProcessing::save_probes(*sim, partial_exporter,true);
    }

    if constexpr (AutoGenerated::FlagCompileTime::use_cereal_serde)
    {
      //Save work with gpu but not loading 
      // do_serde(case_data);
    }

    case_data.simulation->mc_unit.reset();
  }

  std::optional<Core::CaseData> load(const ExecInfo& exec,
                                     const UserControlParameters&& params,
                                     std::optional<Simulation::Feed::SimulationFeed> feed)
  {
    return impl_load(exec, std::forward<decltype(params)>(params), std::move(feed));
  }

} // namespace Core
