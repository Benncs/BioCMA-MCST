#ifndef __GENERATED_MODEL_WRAP__HPP__
#define __GENERATED_MODEL_WRAP__HPP__

#include "common/kokkos_vector.hpp"
#include "mc/domain.hpp"
#include <common/execinfo.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <span>
#include <string_view>

/* ---------------------------------------------------------
 *                  BIOCMA_MC ST
 *
 * This file has been autogenerated by Meson during project
 * configuration.
 *
 * Be careful when modifying values.
 * ---------------------------------------------------------
 */

namespace AutoGenerated
{
  std::unique_ptr<MC::MonteCarloUnit>
  wrap_init_model_selector(int i_model,
                           const ExecInfo& info,
                           size_t number_particle,
                           std::span<double> liq_volume,
                           const NeighborsView<ComputeSpace>& liquid_neighbors,
                           double x0,
                           double& total_mass) noexcept;

  int get_model_index_from_name(std::string_view model_name) noexcept;

  std::vector<std::string> get_model_list() noexcept;
} // namespace AutoGenerated

#endif //__GENERATED_MODEL_WRAP__HPP__
