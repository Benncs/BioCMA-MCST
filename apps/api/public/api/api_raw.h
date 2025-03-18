#ifndef __BIOMC_API_RAW_H__
#define __BIOMC_API_RAW_H__

#include <stddef.h> //NOLINT
#include <stdint.h> //NOLINT
#ifdef __cplusplus
extern "C"
{

  // Forward declaration
  namespace Api
  {
    struct SimulationInstance;
  } // namespace Api
#endif

#ifdef __cplusplus
  typedef struct Api::SimulationInstance* Handle; // NOLINT
#else
typedef struct Opaque*
    Handle; // NOLINT //In C we only need ptr type so Opaque doesn´t need to exist
#endif

  // /**
  //  * @brief Initialize a simulation instance handle for shared memory.
  //  *
  //  * This function creates a simulation instance with a configuration, suitable
  //  * for single-node or single-process execution.
  //  *
  //  * @param id A unique identifier for the simulation instance.
  //  * @param thread_per_process The number of threads allocated per process.
  //  * @return A `Handle` to the simulation instance, or `NULL` if initialization failed.
  //  */
  // Handle init_handle_shared(uint64_t id, uint32_t thread_per_process);

  // // TODO: Enable when USE_MPI is defined
  // /**
  //  * @brief Initialize a raw simulation instance handle with MPI support.
  //  *
  //  * This function creates a simulation instance that is aware of MPI configurations,
  //  * suitable for distributed simulations.
  //  *
  //  * @param n_rank The total number of ranks in the MPI group.
  //  * @param current_rank The rank ID for this instance.
  //  * @param id A unique identifier for the simulation instance.
  //  * @param thread_per_process The number of threads allocated per process.
  //  * @return A `Handle` to the simulation instance, or `NULL` if initialization failed.
  //  */
  // Handle init_handle_raw(int n_rank, int current_rank, uint64_t id, uint32_t thread_per_process);

  Handle init_handle_raw(int argc, char** argv);

  /**
   * @brief Delete a simulation instance handle.
   *
   * This function frees the resources associated with the simulation instance.
   * After calling this function, the `Handle` is no longer valid.
   *
   * @param handle The handle to the simulation instance to delete.
   */
  void delete_handle(Handle* handle);

  void finalize();

  int exec(Handle handle);
  int register_initial_condition(Handle* handle);
  int apply(Handle handle, int to_load);

  /* REGISTER */
  int register_result_path(Handle handle, const char* c);
  int register_cma_path_recursive(Handle handle, const char* c);
  int register_cma_path(Handle handle, const char* c);
  int register_serde(Handle handle, const char* c);
  int register_model_name(Handle handle, const char* c);

  int register_initializer_path(Handle handle, const char* c);
  /*Parameters*/
  // NOLINTBEGIN

  typedef struct wrap_c_param_t
  {
    double biomass_initial_concentration; ///< Initial concentration of biomass.
    double final_time;                    ///< Final time for the simulation (in seconds).
    double delta_time;                    ///< Time step for the simulation (in seconds).
    uint64_t number_particle;             ///< Number of particles in the simulation.
    int32_t n_thread;                     ///< Number of threads to use for simulation.
    uint32_t number_exported_result;      ///< Number of results to be exported.
    int recursive;                        ///< Flag to enable recursive processing.
    int force_override;                   ///< Flag to allow overwriting of existing results.
    int serde;
  } Param;
  // NOLINTEND
  void show_user_param(const Param* params);

  void repr_user_param(const Param* params, char** repr);
  Param make_params(double biomass_initial_concentration,
                    double final_time,
                    double delta_time,
                    uint64_t number_particle,
                    uint32_t number_exported_result);

  int register_parameters(Handle handle, Param* params);

  /*Feed*/
  int set_feed_constant(Handle,
                        double _f,
                        size_t n_species,
                        double* _target,
                        size_t* _species,
                        size_t n_position,
                        size_t* _position,
                        int gas,
                        int fed_batch);
  // /bool set_feed_constant(double _f, std::span<double> _target, std::span<std::size_t> _position,
  // std::span<std::size_t> _species,bool gas=false);

#ifdef __cplusplus
}
#endif

#endif //!__BIOMC_API_RAW_H__
