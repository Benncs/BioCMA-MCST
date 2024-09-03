#include <mc/mcinit.hpp>
#include <string_view>
#include <unordered_map>

@INCLUDES@

    namespace AutoGenerated
{
  std::unique_ptr<MC::MonteCarloUnit> wrap_init_model_selector(
    const int i_model,
      const ExecInfo &info,
      size_t numper_particle,
      std::span<double> liq_volume,
      CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors,
      double x0) noexcept
  {
    switch (i_model)
    {
      @SWITCH_BODY@

          default:
      {
        return MC::init<DefaultModel>(
            info, numper_particle, liq_volume, liquid_neighbors, x0);
      }
    }
  }

  int get_model_index_from_name(std::string_view model_name) noexcept
  {
    @MODEL_INDEX_MAP@
    

    auto it = map.find(model_name);
    if (it != map.end()) {
        return it->second;
    }
    std::cerr<<"Model not found, using DefaultModel instead"<<std::endl;
    return -1;
  }

  std::vector<std::string> get_model_list() noexcept
  {
    return {};
  }
} // namespace AutoGenerated
