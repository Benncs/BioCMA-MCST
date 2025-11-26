#include <algorithm>
#include <cma_utils/alias.hpp>
#include <limits>
namespace CmaUtils
{

  //   double get_min_residence_time(const TransitionnerPtrType& iterator)
  //   {
  //     const std::size_t n_states = iterator->size();
  //     double min_residence_time = std::numeric_limits<double>::max();
  //     for (std::size_t i_state = 0; i_state < n_states; ++i_state)
  //     {
  //       // Naive min-research, as COO is not sorted we can´t do better
  //       without
  //       // convert COO to CSR/CSC or sort
  //       // + function exececuted only once
  //       const auto state = iterator->get_at(i_state);
  //       const auto liquid = state->get_liquid();
  //       const auto coo_matrix = liquid->transition();
  //       const auto liquid_volumes = liquid->volume();
  //       const auto data = coo_matrix->values();
  //       const auto rows = coo_matrix->row_indices();
  //       const auto cols = coo_matrix->col_indices();

  //       for (std::size_t k = 0; k < data.size(); ++k)
  //       {
  //         // Only way to find diagonal with not sorted COO
  //         if (rows[k] == cols[k])
  //         {
  //           const auto volume = liquid_volumes[rows[k]];
  //           if (volume != 0.)
  //           {
  //             double residence_time =
  //                 -data[k] / volume; // minus cause transition has negative
  //                                    // diagonal (leaving flow)
  //             min_residence_time = std::min(residence_time,
  //             min_residence_time);
  //           }
  //         }
  //       }
  //     }
  //     return min_residence_time;
  //   }

  // Does this function should be moved into RCMTool ?
  double get_min_residence_time(const TransitionnerPtrType& iterator) noexcept
  {
    /*
    This new implementation only determine minimum residence time defined as
    min(tau)=min(qi/vi) where qi=sum(flow) in compartment i First impl
    determined the smaller flow min(flowi/vi) which leads to smaller value
    */
    const std::size_t n_states = iterator->size();
    double min_residence_time = std::numeric_limits<double>::max();
    for (std::size_t i_state = 0; i_state < n_states; ++i_state)
    {
      const auto state = iterator->get_at(i_state);
      const auto liquid = state->get_liquid();
      const auto out_flows = liquid->out_flows();
      const auto liquid_volumes = liquid->volume();
      KOKKOS_ASSERT(liquid_volumes.size() == out_flows.size());
      for (std::size_t k = 0; k < out_flows.size(); ++k)
      {
        const auto volume = liquid_volumes[k];
        if (volume > 0.)
        {
          const double residence_time = out_flows[k] / volume;
          min_residence_time = std::min(residence_time, min_residence_time);
        }
      }
    }
    return min_residence_time;
  }

} // namespace CmaUtils