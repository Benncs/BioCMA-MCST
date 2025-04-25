
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <iostream>
#include <mc/prng/prng.hpp>
#include <mc/prng/prng_extension.hpp>

#ifdef CI_TEST
constexpr int N = 4000;
constexpr double tolerance = 0.2;
#else
constexpr int N = 40'000'000;
constexpr double tolerance = 0.05;
#endif

constexpr int test_seed = 1407;
// Distribution parameters
constexpr double mu = 0.0;
constexpr double sigma = 1.0;


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
  MC::KPRNG::pool_type rand_pool(test_seed);
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
  double relative_error{};
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

template <int N, typename T, typename Device, MC::Distributions::ProbabilityLaw<T, Device> F>
void moment_test(F dist, std::string name, T tol = tolerance, int seed = test_seed)
{

  Kokkos::View<T[N]> results("results");
  MC::KPRNG::pool_type rand_pool(seed);

  T sum = 0;
  T sum_sq = 0;
  T sum_skew = 0;
  Kokkos::parallel_reduce(
      name,
      N,
      KOKKOS_LAMBDA(const int i, T& _s, T& _sq, T& _sk) {
        auto rand_gen = rand_pool.get_state();
        const auto val = dist.draw(rand_gen);
        KOKKOS_ASSERT(!Kokkos::isnan(val) && !Kokkos::isinf(val));
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

  T mean = sum / static_cast<T>(N);
  T var = sum_sq / static_cast<T>(N) - mean * mean;
  T skew = (sum_skew / N - 3 * mean * var - mean * mean * mean) / std::pow(var, 1.5);

  check_moment(mean, dist.mean(), tol, "Mean");
  check_moment(var, dist.var(), tol, "Variance");
  check_moment(skew, dist.skewness(), tol, "Skewness");
}

template <int N> void test_norminv()
{
  // Std normal
  double mu = 0.;
  double sigma = 1;
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
  MC::KPRNG::pool_type rand_pool(test_seed);
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
  test_norminv<N>();
  std::cerr << "test with " << N << " samples" << std::endl;
  std::cerr << "Normal" << std::endl;
  moment_test<N, double, Device>(MC::Distributions::Normal<double>{0, 1}, "Normal");

  std::cerr << "Normal2" << std::endl;

  moment_test<N, double, Device>(MC::Distributions::Normal<double>{mu, sigma}, "Normal2");

  MC::Distributions::LogNormal<double> log(mu, sigma);
  std::cerr << "Log-Normal" << std::endl;
  moment_test<N, double, Device>(log, "Log-Normal");
  std::cerr << "Skew-Normal" << std::endl;
  double xi = 0;
  double omega = 1.;
  double alpha = 0;
  moment_test<N, double, Device>(MC::Distributions::SkewNormal<double>{xi, omega, alpha},
                                 "Skew-Normal");

  std::cerr << "Truncated-Normal" << std::endl;

  moment_test<N, double, Device>(MC::Distributions::TruncatedNormal<double>{1., 0.33, 0., 5.},
                                 "Truncated-Normal");

  std::cerr << "Truncated-Normal2" << std::endl;
  moment_test<N, double, Device>(MC::Distributions::TruncatedNormal<double>(-5., 1.4, -10., 1.),
                                 "Truncated-Normal2");

  std::cerr << "Truncated-Normal3" << std::endl;
  for (int i = 0; i < 5; ++i)
  {
    double factor = 10000.;
    double mu_t = 4e-5 * factor;
    double maxsample = 9e-5 * factor;
    moment_test<N, double, Device>(
        MC::Distributions::TruncatedNormal<double>(mu_t, 0.01, 0., maxsample),
        "Truncated-Normal3",
        0.3,
        i * test_seed);
  }

  std::cerr << "Scaled double Truncated-Normal3" << std::endl;

  constexpr float l_min_um = 0.9;
  constexpr float l_max_um = 5;
  constexpr auto length_c_dist =
      MC::Distributions::TruncatedNormal<double>(l_min_um, l_min_um / 5., 0.5 * l_min_um, l_max_um);

  for (int i = 0; i < 5; ++i)
  {
    moment_test<N, double, Device>(
        length_c_dist, "double Truncated-Normal3", 0.2, AutoGenerated::debug_MC_RAND_DEFAULT_SEED);
  }

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
