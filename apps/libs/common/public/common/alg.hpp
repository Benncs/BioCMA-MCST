#ifndef __COMMON_ALG_HPP__
#define __COMMON_ALG_HPP__
#include <Kokkos_Sort.hpp>
#include <common/common.hpp>
#include <functional>
#include <limits>
#include <sorting/Kokkos_BinSortPublicAPI.hpp>
#include <sorting/impl/Kokkos_CopyOpsForBinSortImpl.hpp>
#include <sorting/impl/Kokkos_SortByKeyImpl.hpp>
double naive_newton(const std::function<double(double)>& f,
                    double x_init,
                    bool* success,
                    double tolerance = std::numeric_limits<double>::epsilon());

namespace Sorting
{
  namespace
  {
    template <typename PermType, typename... Views>
    void
    copy_data(std::size_t n_allocated_elements,
              std::size_t n_used_elements,
              const PermType& perm,
              std::size_t src_offset,
              Views... views)
    {
      // The following may not work on GPU

      // auto tmp_views = std::make_tuple(Kokkos::create_mirror(views)...);

      // Kokkos::parallel_for(
      //     "Sort::copy_data",
      //     Kokkos::RangePolicy<>(0, n_allocated_elements),
      //     KOKKOS_LAMBDA(int i) {
      //       if (perm(i) >= n_used_elements)
      //       {
      //         return;
      //       }

      //       const int idx = src_offset + perm(i);

      //       std::apply(
      //           [&](auto&... tmp)
      //           {
      //             ((Kokkos::Impl::CopyOp<
      //                  std::decay_t<decltype(views)>,
      //                  std::decay_t<decltype(views)>>::copy(tmp, i, views,
      //                  idx)),
      //              ...);
      //           },
      //           tmp_views);
      //     });

      // auto do_deep_copy = [](auto& dst, auto& src)
      // { Kokkos::deep_copy(dst, src); };

      // std::apply([&](auto&... tmps) { (do_deep_copy(views, tmps), ...); },
      //            tmp_views);
    }

  } // namespace

  template <typename RefView, typename... ViewType>
  void
  sort_soa(std::size_t n0,
           std::size_t nn,
           std::size_t n_e,
           std::size_t n_max,
           RefView ref_sorting,
           ViewType... views)
  {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    PROFILE_SECTION("SORT")
    if (nn - n0 == 1 || nn < n0)
    {
      return;
    }
    const int bin_size = 2048;
    using view_type = decltype(ref_sorting);
    auto binop = Kokkos::BinOp1D<view_type>(bin_size, n0, nn);

    auto sorter = Kokkos::BinSort<view_type, decltype(binop)>(
        ExecSpace(), ref_sorting, 0, n_e, binop, false);

    sorter.create_permute_vector();
    auto perm = sorter.get_permute_vector();
    auto src_offset = 0;
    copy_data(n_max, n_e, perm, src_offset, views...);
  }
} // namespace Sorting

#endif //__COMMON_ALG_HPP__