#include <common/common.hpp>

SimulationParameters SimulationParameters::m_default()
{
  SimulationParameters params;
  params.n_species = 1;
  params.flow_files = {"./cma_data/"};
  params.verbose = false;
  params.n_compartments = 0;
  params.t_per_flow_map = 0;
  params.user_params = UserControlParameters::m_default();

  return params;
}

UserControlParameters UserControlParameters::m_default()
{
  return {"", 
  false, 0, 
  0, 0, "", -1, 0,"model_light"};
}
