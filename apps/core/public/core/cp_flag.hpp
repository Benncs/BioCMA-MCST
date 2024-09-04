#ifndef __CORE_COMPILE_TIME_FLAGS_HPP__
#define __CORE_COMPILE_TIME_FLAGS_HPP__

namespace FlagCompileTIme
{
#ifndef __FLAG_APP_VERBOSE__
  constexpr bool verbose = false;
#else
  constexpr bool verbose = true;
#endif

#ifdef NO_MPI
  constexpr bool use_mpi = false;
#else
  constexpr bool use_mpi = true;
#endif
} // namespace FlagCompileTIme

#endif
