#ifndef __PMESSAGE_HPP__
#define __PMESSAGE_HPP__


#include <messages/mpi_types.hpp>


namespace MPI_W
{
  static bool is_mpi_init = false;
  void critical_error() noexcept;
  void barrier()noexcept;
  void finalize()noexcept;
  enum class SIGNALS
  {
    STOP,
    RUN
  };

} // namespace MPI_W

#endif //__PMESSAGE_HPP__