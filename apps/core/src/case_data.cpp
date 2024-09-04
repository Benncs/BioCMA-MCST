#include <core/case_data.hpp>
#include <host_specific.hpp>
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

    auto *const sim = case_data.simulation.release();

    f_run(case_data.exec_info,
          std::move(*sim),
          case_data.params,
          std::move(case_data.transitioner));
  }

} // namespace Core