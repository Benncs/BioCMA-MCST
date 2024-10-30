#ifndef __EXEC_INFO_HPP__
#define __EXEC_INFO_HPP__

#include <cstddef>
#include <iostream>
#include <cstdint> 

struct ExecInfo
{
  uint32_t n_rank;
  uint32_t current_rank;
  uint32_t thread_per_process;
  bool verbose;
  uint64_t run_id;
#if defined(__cpp_lib_hardware_interference_size)
  // default cacheline size from runtime
  static constexpr size_t cache_line_size =
      std::hardware_destructive_interference_size;
#else
  // most common cacheline size otherwise
  static constexpr size_t cache_line_size = 64;
#endif
};



inline std::ostream &operator<<(std::ostream &stream, const ExecInfo &obj)
{
  stream << obj.run_id << "\t" << obj.n_rank << "\t" << obj.thread_per_process
         << "\t";
  return stream;
}

#endif //__EXEC_INFO_HPP__