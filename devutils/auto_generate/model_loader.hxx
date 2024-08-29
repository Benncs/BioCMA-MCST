#ifndef __GENERATED_MODEL_WRAP__HPP__
#define __GENERATED_MODEL_WRAP__HPP__

#include <common/execinfo.hpp>
#include <mc/mcinit.hpp>
#include <memory>
#include <span>
#include <string_view>

namespace AutoGenerated
{
  std::unique_ptr<MC::MonteCarloUnit> wrap_init_model_selector(
      const ExecInfo &info,
      size_t numper_particle,
      std::span<double> liq_volume,
      CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors,
      double x0);

  int get_model_index_from_name(std::string_view model_name);


  std::vector<std::string> get_model_list();
} // namespace AutoGenerated

#endif //__GENERATED_MODEL_WRAP__HPP__
