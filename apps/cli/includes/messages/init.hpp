#ifndef __APP_INIT_HPP__
#define __APP_INIT_HPP__

#include <common/execinfo.hpp>
#include <common/execinfo.hpp>
#include <mpi.h>
#include <messages/pmessage.hpp>
#include <span> 

ExecInfo init_mpi(int argc, char **argv);


template <typename... Args>
void host_dispatch(const ExecInfo &info, MPI_SIGNALS &&sign, Args &&...args)
{

  // Send message to all other processes
  for (int j = 1; j < static_cast<int>(info.n_rank); ++j)
  {
    MPI_Send(&sign, sizeof(sign), MPI_CHAR, j, 0, MPI_COMM_WORLD);
    (
        [&]<typename T>(T &&arg)
        {
          size_t s;
          void *buf;
          if constexpr (std::is_same_v<std::decay_t<T>, std::span<double>>)
          {
            s = arg.size();
            buf = arg.data();
            MPI_Send(&s, 1, MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD);
          }
          else
          {
            s = sizeof(T);
            buf = &arg;
          }

          MPI_Send(buf, s, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
        }(std::forward<Args>(args)),
        ...);
  }
}

#endif //__APP_INIT_HPP__