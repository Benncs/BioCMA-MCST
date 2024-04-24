#ifndef __EXEC_INFO_HPP__
#define __EXEC_INFO_HPP__

#include <cstddef>
#include <iostream>

struct ExecInfo
{
  size_t n_rank;
  size_t current_rank;
  size_t thread_per_process;
  bool verbose;
  size_t run_id;
};

inline std::ostream &operator<< (std::ostream &stream, const ExecInfo & obj) {
    stream<<obj.run_id<<"\t"<<obj.n_rank<<"\t"<<obj.thread_per_process<<"\t";
    return stream;
}

#endif //__EXEC_INFO_HPP__