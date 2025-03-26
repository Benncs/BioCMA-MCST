#include "Kokkos_Core.hpp"
#include <biocma_cst_config.hpp>
#include <common/execinfo.hpp>
#include <core/case_data.hpp>
#include <core/global_initaliser.hpp>
#include <filesystem>
#include <host_specific.hpp>
#include <impl_post_process.hpp>
#include <optional>
#include <serde.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <worker_specific.hpp>

#ifdef NO_MPI
namespace WrapMPI
{
  inline bool is_initialized() noexcept
  {
    return false;
  }
} // namespace WrapMPI
#else
#  include "mpi_w/message_t.hpp"
#  include <mpi_w/wrap_mpi.hpp>
#endif
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
    //TODO Find orphan profile_section : Kokkos::Profiling::popRegion() called outside  of an actve region. Previous regions: 
    
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
      PostProcessing::save_probes(*sim, partial_exporter, true);
    }

    if constexpr (AutoGenerated::FlagCompileTime::use_cereal_serde)
    {
      // Save work with gpu but not loading
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

  ExecInfo runtime_init(int argc, char** argv, std::optional<std::size_t> force_run_id)
  {
    ExecInfo info{};

    int rank = 0;
    int size = 1;

    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
#ifndef NO_MPI
      if (!WrapMPI::is_initialized())
      {
        MPI_Init(&argc, &argv);
      }
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    }

    if (!Kokkos::is_finalized() && !Kokkos::is_initialized())
    {
      Kokkos::initialize(argc, argv);
    }

#ifndef NO_MPI
    WrapMPI::barrier();
#endif
    if (rank == 0)
    {
      Kokkos::print_configuration(std::cout);
      if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
      {
        std::cout << "MPI NODES: " << size << std::endl;
      }
    }

    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
#ifndef NO_MPI
      std::atexit(WrapMPI::finalize);
#endif
    }

    info.current_rank = static_cast<size_t>(rank);
    info.n_rank = static_cast<size_t>(size);

    info.thread_per_process = Kokkos::num_threads();

    const auto id_seed =
        (force_run_id.has_value())
            ? *force_run_id
            : static_cast<size_t>(time(nullptr) * info.n_rank * info.thread_per_process);

    std::string s = std::to_string(id_seed);
    info.run_id = std::hash<std::string>{}(s);

    return info;
  }

} // namespace Core
