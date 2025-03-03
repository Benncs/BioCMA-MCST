#include "Kokkos_Assert.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_Random.hpp"
#include "decl/Kokkos_Declare_OPENMP.hpp"
#include <Kokkos_Core.hpp>
#include <cassert>
#include <iostream>
#include <limits>
#include <mc/prng/prng_extension.hpp>
#include <utility>

constexpr int test_seed = 1407;
// Distribution parameters
constexpr double mu = 0.0;
constexpr double sigma = 1.0;

constexpr double tolerance = 0.05;

template <typename T> void basic_test()
{
  T xi = 0;
  T omega = 1.;
  T alpha = 0;
  const int N = 1000;
  Kokkos::View<T*> log_results("log_results", N);
  Kokkos::View<T*> skew_results("skew_results", N);
  MC::Distributions::LogNormal<T> log(mu, sigma);
  MC::Distributions::SkewNormal<T> skewnormal{xi, omega, alpha};
  Kokkos::Random_XorShift64_Pool<> rand_pool(test_seed);
  Kokkos::parallel_for(
      "TestKernel", N, KOKKOS_LAMBDA(const int i) {
        auto rand_gen = rand_pool.get_state();

        skew_results(i) = skewnormal.draw(rand_gen);
        log_results(i) = log.draw(rand_gen);
        rand_pool.free_state(rand_gen);
      });
  Kokkos::fence();
}
void check_moment(double empirical, double theoretical, double tol, const char* name)
{
  double relative_error;
  if (theoretical != 0)
  {
    relative_error = std::abs(empirical - theoretical) / theoretical;
  }
  else
  {
    relative_error = std::abs(empirical);
  }
  std::cerr << "Test  " << name << "\n";
  std::cerr << "  Empirical: " << empirical << "\n";
  std::cerr << "  Theoretical: " << theoretical << "\n";
  std::cerr << "  Relative error: " << relative_error << "\n";
  if (std::abs(relative_error) >= tol)
  {
    std::cerr << "Failed  " << name << "\n";
    assert(false);
  }
}

template <int N, typename Device, MC::Distributions::ProbabilityLaw<double, Device> F>
void moment_test(F dist, std::string name)
{

  Kokkos::View<double[N]> results("results");
  Kokkos::Random_XorShift64_Pool<> rand_pool(test_seed);

  double sum = 0;
  double sum_sq = 0;
  double sum_skew = 0;
  Kokkos::parallel_reduce(
      name,
      N,
      KOKKOS_LAMBDA(const int i, double& _s, double& _sq, double& _sk) {
        auto rand_gen = rand_pool.get_state();
        const auto val = dist.draw(rand_gen);
        results(i) = val;
        _s += val;
        _sq += (val * val);
        _sk += (val * val * val);
        rand_pool.free_state(rand_gen);
      },
      sum,
      sum_sq,
      sum_skew);
  Kokkos::fence();

  double mean = sum / static_cast<double>(N);
  double var = sum_sq / static_cast<double>(N) - mean * mean;
  double skew = (sum_skew / N - 3 * mean * var - mean * mean * mean) / std::pow(var, 1.5);

  check_moment(mean, dist.mean(), tolerance, "Mean");
  check_moment(var, dist.var(), tolerance, "Variance");
  check_moment(skew, dist.skewness(), tolerance, "Skewness");
}

void test_norminv()
{
  // Std normal
  double mu = 0.;
  double sigma = 1;
  const int N = 10000;
  // double x_min=-10;
  // double x_max = 10;
  // double dx = (x_max-x_min)/N;
  // Kokkos::View<double*> results2("results", N);
  Kokkos::View<double*> results("results", N);
  Kokkos::Random_XorShift64_Pool<> rand_pool(test_seed);

  Kokkos::parallel_for(
      "Generate and test norminv", N, KOKKOS_LAMBDA(const int i) {
        auto rand_gen = rand_pool.get_state();
        double u = rand_gen.drand(0.0, 1.0);
        results(i) = MC::Distributions::norminv(u, mu, sigma);

        // double x = x_min+i*dx;
        // results2(i) = MC::Distributions::norminv(x, mu, sigma);

        rand_pool.free_state(rand_gen);
      });

  Kokkos::fence();

  auto host_results = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(host_results, results);

  double mean = 0.0;
  double variance = 0.0;
  for (int i = 0; i < N; i++)
  {
    mean += host_results(i);
  }
  mean /= N;

  for (int i = 0; i < N; i++)
  {
    variance += (host_results(i) - mean) * (host_results(i) - mean);
  }
  variance /= (N - 1);

  assert(std::abs(mean) < 0.1);
  assert(std::abs(variance - 1.0) < 0.1);

  
  // Kokkos::parallel_for(
  //     "test norminv", Kokkos::RangePolicy(1,N-1), KOKKOS_LAMBDA(const int i) {
  //       Kokkos::printf("%lf\r\n",results2(i));
  //       KOKKOS_ASSERT(results2(i-1)<=results2(i));
  //     });

}

template <int N> void test_truncated()
{
  MC::Distributions::TruncatedNormal<double> dist(2, 2. / 3., 0., 4.);
  Kokkos::Random_XorShift64_Pool<> rand_pool(test_seed);
  Kokkos::parallel_for(
      "test_truncated", N, KOKKOS_LAMBDA(const int _) {
        auto gen = rand_pool.get_state();
        auto val = dist.draw(gen);
        KOKKOS_ASSERT(val >= 0.); // bounds should be ]a,b[ but test in [a,b]
        KOKKOS_ASSERT(val <= 4.)
        rand_pool.free_state(gen);
      });
}

int main()
{
  using Device = Kokkos::DefaultExecutionSpace;
  Kokkos::initialize();
  Kokkos::print_configuration(std::cout);
  basic_test<float>();
  basic_test<double>();
  test_norminv();

  constexpr int N = 40'000'000;
  std::cerr << "Normal" << std::endl;
  moment_test<N, Device>(MC::Distributions::Normal<double>{0, 1}, "Normal");

  std::cerr << "Normal2" << std::endl;

  moment_test<N, Device>(MC::Distributions::Normal<double>{mu, sigma}, "Normal2");

  MC::Distributions::LogNormal<double> log(mu, sigma);
  std::cerr << "Log-Normal" << std::endl;
  moment_test<N, Device>(log, "Log-Normal");
  std::cerr << "Skew-Normal" << std::endl;
  double xi = 0;
  double omega = 1.;
  double alpha = 0;
  moment_test<N, Device>(MC::Distributions::SkewNormal<double>{xi, omega, alpha}, "Skew-Normal");

  std::cerr << "Truncated-Normal" << std::endl;

  moment_test<N, Device>(MC::Distributions::TruncatedNormal<double>{1., 0.33, 0., 5.},
                         "Truncated-Normal");

  std::cerr << "Truncated-Normal2" << std::endl;
  moment_test<N, Device>(MC::Distributions::TruncatedNormal<double>(-5., 1.4, -10., 1.),
                         "Truncated-Normal2");
  test_truncated<N>();
  std::cerr << "Skew-Normal2" << std::endl;

  // Dont work
  //  xi = 3.;
  //  omega = 3./7.;
  //  alpha = 4;
  //  moment_test<N,Device>(MC::SkewNormal<double>{xi, omega, alpha},"Skew-Normal2");
  //  std::cerr << "Skew-Normal3" << std::endl;
  //  xi = -1.;
  //  omega = 1.;
  //  alpha = -7;
  //  moment_test<N,Device>(MC::SkewNormal<double>{xi, omega, alpha},"Skew-Normal3");

  Kokkos::finalize();
  return 0;
}
