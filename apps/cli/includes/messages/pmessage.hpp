#ifndef __PMESSAGE_HPP__
#define __PMESSAGE_HPP__

#include <concepts>
#include <cstddef>


struct InitMessage
{
  size_t n_compartments;
  size_t n_neighbor;
  double d_t;
};

enum class MPI_SIGNALS
{
  STOP,
  RUN
};




#endif //__PMESSAGE_HPP__