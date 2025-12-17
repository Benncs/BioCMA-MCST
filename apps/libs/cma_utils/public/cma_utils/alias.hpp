#ifndef __CMA_UTILS_ALIAS_HPP__
#define __CMA_UTILS_ALIAS_HPP__

#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <lib.rs.h>

namespace CmaUtils
{
  /** @brief Opaque type for flowmap transitioner
    @note: underlying type is rust box which has to be considered as a
    non-movable unique_ptr
  */
  using TransitionnerPtrType = rust::Box<TransitionerWrapper>;

  /** @brief Opaque type for iteration state
    @note: underlying type is rust box which has to be considered as a
    non-movable unique_ptr
  */
  using IterationStatePtrType = ::rust::Box< ::IterationStateWrapper>;

  /** @brief Opaque type for Naive COO matrix
    @note: underlying type is rust box which has to be considered as a
    non-movable unique_ptr
  */
  using StateCooMatrixType = ::rust::Box< ::CooMatrixWrap>;

  /** @brief Determine the smallest compartment residence time
      For cases with multiple flowmaps, all the different states are taken into
     account*/
  double get_min_residence_time(const TransitionnerPtrType& iterator) noexcept;

} // namespace CmaUtils

#endif
