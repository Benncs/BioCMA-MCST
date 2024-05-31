#include <common/common.hpp>

SimulationParameters SimulationParameters::m_default()
{
  SimulationParameters params;
  params.n_particles = 0;
  params.n_species = 1;
  params.final_time = 0;
  params.flow_files = {"./cma_data/"};
  params.d_t = 0;
  params.n_threads = -1;
  params.verbose = false;
  params.n_compartments = 0;
  return params;
}
