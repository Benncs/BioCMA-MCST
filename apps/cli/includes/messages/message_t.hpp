#ifndef __PMESSAGE_HPP__
#define __PMESSAGE_HPP__


namespace MPI_W
{
  // static bool is_mpi_init = false;
  int critical_error() noexcept;
  void barrier() noexcept;
  void finalize() noexcept;
  enum class SIGNALS
  {
    STOP,
    RUN,
    NOP
  };



} // namespace MPI_W

#endif //__PMESSAGE_HPP__