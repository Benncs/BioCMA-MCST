#include "Kokkos_Macros.hpp"
#include "Kokkos_MathematicalFunctions.hpp"
#include "common/env_var.hpp"
#include <models/config_loader.hpp>

#include <models/fixed_length.hpp>

#include <Kokkos_sampling/metropolis.hpp>

namespace Models
{
  FixedLength::Config
  FixedLength::get_config(const ExecInfo& info, const std::size_t n)
  {
    using float_type = Self::FloatType;
    Kokkos::View<float_type*, ComputeSpace> samples("samples", n);

    float_type lambda = Kokkos::log(2.F) / 1e-6F;

    lambda = Common::read_env_or("VLAMBDA", lambda);

    const std::size_t nbin = Common::read_env_or("VNBIN", 50);

    if (lambda != 0.F)
    {
      const size_t n_rank = info.n_rank;

      const size_t current_rank = info.current_rank;

      std::size_t Ntot
          = Common::read_env_or("__N_TOTAL_PARTICLE__", n * n_rank);

      const std::size_t n_edges = nbin + 1;
      const float_type dl
          = (Self::l_max_m - Self::l_min_m) / static_cast<float_type>(nbin);

      std::vector<double> e_edge_h(n_edges);
      for (std::size_t i = 0; i < n_edges; ++i)
      {
        const double l_i = static_cast<double>(Self::l_min_m)
                           + static_cast<double>(i) * static_cast<double>(dl);
        e_edge_h[i] = std::exp(-static_cast<double>(lambda) * l_i);
      }
      const double e_range = e_edge_h[0] - e_edge_h[n_edges - 1];

      std::vector<double> n_real(nbin);
      std::vector<std::size_t> n_floor(nbin);
      double sum_floor = 0.0;

      for (std::size_t i = 0; i < nbin; ++i)
      {
        const double mass_i = (e_edge_h[i] - e_edge_h[i + 1]) / e_range;
        n_real[i] = mass_i * static_cast<double>(Ntot);
        n_floor[i] = static_cast<std::size_t>(std::floor(n_real[i]));
        sum_floor += static_cast<double>(n_floor[i]);
      }

      const auto remainder = Ntot - static_cast<std::size_t>(sum_floor);

      std::vector<std::size_t> order(nbin);
      std::iota(order.begin(), order.end(), 0);
      std::sort(
          order.begin(),
          order.end(),
          [&](std::size_t a, std::size_t b)
          { return (n_real[a] - n_floor[a]) > (n_real[b] - n_floor[b]); });

      std::vector<std::size_t> n_bin_global_h = n_floor;
      for (std::size_t k = 0; k < remainder; ++k)
      {
        n_bin_global_h[order[k]] += 1;
      }

      std::vector<std::size_t> bin_start_global_h(nbin + 1);
      bin_start_global_h[0] = 0;
      for (std::size_t i = 0; i < nbin; ++i)
      {
        bin_start_global_h[i + 1] = bin_start_global_h[i] + n_bin_global_h[i];
      }

      Kokkos::View<double*, ComputeSpace> e_edge("e_edge", n_edges);
      Kokkos::View<const double*, Kokkos::HostSpace> e_edge_hv(e_edge_h.data(),
                                                               n_edges);
      Kokkos::deep_copy(e_edge, e_edge_hv);

      Kokkos::View<std::size_t*, Kokkos::HostSpace> bin_start_host(
          bin_start_global_h.data(), nbin + 1);
      Kokkos::View<std::size_t*, ComputeSpace> bin_start_global
          = Kokkos::create_mirror_view_and_copy<ComputeSpace>(ComputeSpace(),
                                                              bin_start_host);

      Kokkos::RangePolicy<ComputeSpace> pl0(0, n);
      Kokkos::parallel_for(
          "init_fixed_length_exact", pl0, KOKKOS_LAMBDA(const std::size_t j) {
            // const std::size_t k = rank_global_offset + j;
            //
            const std::size_t k = current_rank + j * n_rank;

            // binary search: find largest bin index `bin` such that
            // bin_start_global(bin) <= k < bin_start_global(bin+1)
            std::size_t lo = 0;
            std::size_t hi = nbin; // number of bins
            while (lo + 1 < hi)
            {
              std::size_t mid = (lo + hi) / 2;
              if (bin_start_global(mid) <= k)
              {
                lo = mid;
              }
              else
              {
                hi = mid;
              }
            }
            const std::size_t bin = lo;

            const std::size_t bin_start = bin_start_global(bin);
            const std::size_t bin_stop = bin_start_global(bin + 1);
            const std::size_t n_i_glob = bin_stop - bin_start;
            const std::size_t pos_in_bin = k - bin_start;

            // bin's low-x edge / high-x edge
            const float_type e_hi = static_cast<float_type>(e_edge(bin));
            const float_type e_lo = static_cast<float_type>(e_edge(bin + 1));
            const float_type x_lo
                = Self::l_min_m + static_cast<float_type>(bin) * dl;
            const float_type x_hi
                = Self::l_min_m + static_cast<float_type>(bin + 1) * dl;

            float_type t = (static_cast<float_type>(pos_in_bin) + 0.5F)
                           / static_cast<float_type>(n_i_glob);

            float_type e_val = e_hi - t * (e_hi - e_lo);
            float_type x = -(1.0F / lambda) * Kokkos::log(e_val);
            x = Kokkos::fmax(x_lo, Kokkos::fmin(x_hi, x));
            samples(j) = x;
          });
    }
    else
    {
      throw std::runtime_error("Unimplemented yet");
    }

    return samples;
  }

  // FixedLength::Config
  // FixedLength::get_config(const std::size_t n)
  // {
  //   using float_type = Self::FloatType;
  //   float_type lambda = Kokkos::log(2.F) / 1e-6F;
  //   char* lambda_env = std::getenv("VLAMBDA");

  //   if (lambda_env != nullptr)
  //   {
  //     lambda = static_cast<float_type>(std::stod(lambda_env));
  //     Kokkos::printf("[Config] use env value %f\r\n", lambda);
  //   }

  //   int rc = 0;
  //   Kokkos::View<float_type*, ComputeSpace> samples("samples", n);
  //   if (lambda != 0.)
  //   {
  //     auto target = KOKKOS_LAMBDA(const FixedLength::FloatType x)
  //     {
  //       return lambda * Kokkos::exp(-lambda * x);
  //     };

  //     rc = Sampling::metropolis(target, samples, Self::l_min_m,
  //     Self::l_max_m);
  //   }
  //   else
  //   {
  //     constexpr FixedLength::FloatType mu
  //         = 1.3e-6; // maxlength is 2 and min i 1 so let be between
  //     constexpr FixedLength::FloatType sigma = 0.1;

  //     auto target_distribution = KOKKOS_LAMBDA(const FixedLength::FloatType
  //     x)
  //     {
  //       if (x <= 0.0)
  //       {
  //         return 0.0;
  //       }
  //       auto coefficient = 1.0 / (x * sigma * std::sqrt(2.0 * M_PI));
  //       auto exponent = -std::pow(std::log(x) - mu, 2) / (2.0 * sigma *
  //       sigma);

  //       return coefficient * std::exp(exponent);
  //     };

  //     rc = Sampling::metropolis(
  //         target_distribution, samples, Self::l_min_m, Self::l_max_m);
  //   }

  //   if (rc != 0)
  //   {
  //     throw std::runtime_error("FixedLength init: Error when sampling");
  //   }

  //   return samples;
  // }

} // namespace Models
