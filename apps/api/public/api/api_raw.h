#ifndef __BIOMC_API_RAW_H__
#define __BIOMC_API_RAW_H__

#include <cstddef>
#include <cstdint>
#include <stdint.h> //NOLINT
#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct Handle Handle; // NOLINT

  Handle* init_handle_raw(int n_rank, int current_rank, uint64_t id);
  Handle* load_handle(int n_rank, int current_rank);
  void delete_handle(Handle* handle);
  int exec(Handle* handle);
  int register_initial_condition(Handle* handle);
  void register_parameters(Handle* handle);
  void apply(Handle* handle, int to_load);

  /* REGISTER */
  int register_result_path(Handle* handle, const char* c);
  int register_result_path_recursive(Handle* handle, const char* c);
  int register_cma_path(Handle* handle, const char* c);

#ifdef __cplusplus
}
#endif

#endif //!__BIOMC_API_RAW_H__
