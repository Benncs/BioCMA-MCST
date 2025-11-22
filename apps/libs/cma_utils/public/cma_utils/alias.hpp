#ifndef __CMA_UTILS_ALIAS_HPP__
#define __CMA_UTILS_ALIAS_HPP__

#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <lib.rs.h>

namespace CmaUtils
{
  using TransitionnerPtrType = rust::Box<TransionnerWrapper>;
  using IterationStatePtrType = ::rust::Box<::IterationStateWrapper>;
  using StateCooMatrixType = ::rust::Box<::CooMatrixWrap>;


} // namespace CmaUtils

#endif
