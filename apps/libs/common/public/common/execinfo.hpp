#ifndef __EXEC_INFO_HPP__
#define __EXEC_INFO_HPP__

#include <cstddef>

struct ExecInfo
{
  size_t n_rank;
  size_t current_rank;
  size_t thread_per_process;
  bool verbose;
};

#endif //__EXEC_INFO_HPP__