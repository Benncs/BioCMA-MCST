#ifndef __BIOMC_API_RAW_H__
#define __BIOMC_API_RAW_H__

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct Handle Handle; //NOLINT

  Handle *init_handle_raw(int n_rank,int current_rank);
  Handle* load_handle(int n_rank,int current_rank);
  void delete_handle(Handle *handle);
  int exec(Handle *handle);
  int register_initial_condition(Handle *handle);
  void register_parameters(Handle *handle);
  void apply(Handle *handle);
    // void register_parameters(Handle *handle, Core::SimulationParameters *params);
  // bool register_feed(Handle *handle, Simulation::Feed::SimulationFeed *feed);
#ifdef __cplusplus
}
#endif

#endif //!__BIOMC_API_RAW_H__
