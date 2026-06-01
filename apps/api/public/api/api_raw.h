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

  namespace Simulation::Feed
  {
    struct FeedDescriptor;
  }; // namespace Simulation::Feed

  typedef struct Api::SimulationInstance* Handle;              // NOLINT
  typedef struct Simulation::Feed::FeedDescriptor* FeedHandle; // NOLINT

#else
// In C we only need ptr type so Opaque doesn´t need to exist
typedef struct Opaque* Handle; // NOLINT

// In C we only need ptr type so Opaque doesn´t need to exist
typedef struct OpaqueFeed* FeedHandle; // NOLINT
#endif

  /*FFI Feed descritptor*/
  FeedHandle new_constant_feed_descriptor(double flow, uint64_t input_position);

  int add_species(FeedHandle fh, double value, uint64_t i_species);

  int set_fedbatch(FeedHandle fh);

  int set_output_position(FeedHandle fh, uint64_t output_position);

  int delete_constant_feed_descriptor(FeedHandle* fd);

  int version_is_compatible(int major, int minor, int dev);

  /*FFI Parameters*/

  // NOLINTBEGIN
  typedef struct wrap_c_param_t
  {
    double biomass_initial_concentration; ///< Initial concentration of biomass.
    double final_time;        ///< Final time for the simulation (in seconds).
    double delta_time;        ///< Time step for the simulation (in seconds).
    uint64_t number_particle; ///< Number of particles in the simulation.
    int32_t n_thread;         ///< Number of threads to use for simulation.
    uint32_t number_exported_result; ///< Number of results to be exported.
    int force_override; ///< Flag to allow overwriting of existing results.
    int load_serde;
    int save_serde;
    int uniform_particle_init;
  } Param;
  // NOLINTEND

  void show_user_param(const Param* params);

  void repr_user_param(const Param* params, char** repr);
  Param make_params(double biomass_initial_concentration,
                    double final_time,
                    double delta_time,
                    uint64_t number_particle,
                    uint32_t number_exported_result,
                    int save);

  Param* make_params_ptr(double biomass_initial_concentration,
                         double final_time,
                         double delta_time,
                         uint64_t number_particle,
                         uint32_t number_exported_result,
                         int save);

  void delete_params(Param** params);

  // FFI API

  /**
 * @brief Initialize a simulation instance handle
   *
   * This function creates a simulation instance with a configuration
   *
   * @param argc.
   * @param argv.
   * @return A `Handle` to the simulation instance, or `NULL` if
   initialization failed.
   */
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

  void get_model_list(char** names, const int* n_model);
  void free_model_list(char** names, int n_model);
  // void finalize(); //Do not use it

  int n_rank(Handle);
  int i_rank(Handle);

  int exec(Handle);

  int apply(Handle, int to_load);

  /* REGISTER */
  int register_initial_condition(Handle*);
  int register_result_path(Handle, const char* c);
  int register_cma_path(Handle, const char* c);
  int register_serde(Handle, const char* c);
  int register_model_name(Handle, const char* c);
  int register_parameters(Handle, Param* params);
  int register_initializer_path(Handle, const char* c);

  int set_scalar_buffer(
      Handle, uint64_t rows, uint64_t cols, double* liquid, double* gas_ptr);

  int set_feed_constant(Handle,
                        double flow,
                        double concentration,
                        size_t species,
                        size_t position,
                        int output_position,
                        int gas,
                        int fed_batch);

  int add_feed_descriptor(Handle, FeedHandle, int gas);

#ifdef __cplusplus
}
#endif

#endif //!__BIOMC_API_RAW_H__
