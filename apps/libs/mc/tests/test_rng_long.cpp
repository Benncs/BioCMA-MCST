#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <chrono>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <mc/prng/prng.hpp>
#include <omp.h>
#include <random>
#include <traits/Kokkos_IterationPatternTrait.hpp>
void calculate_correlation(const std::vector<double> &samples, size_t n_sample)
{
  // Test de corrélation
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_xy = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;
  for (size_t i = 0; i < n_sample - 1; ++i)
  {
    double x = samples[i];
    double y = samples[i + 1];

    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_x2 += x * x;
    sum_y2 += y * y;
  }

  double correlation = (n_sample * sum_xy - sum_x * sum_y) /
                       std::sqrt((n_sample * sum_x2 - sum_x * sum_x) *
                                 (n_sample * sum_y2 - sum_y * sum_y));
  std::cout << "Pearson correlation: " << correlation << std::endl;
}

void mean_variance(const std::vector<double> &samples, size_t n_sample)
{
  double mean =
      std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();

  double variance = 0.0;
  for (double sample : samples)
  {
    variance += (sample - mean) * (sample - mean);
  }
  variance /= samples.size();

  std::cout << "Mean : " << mean << std::endl;
  std::cout << "Variance : " << variance << std::endl;
}

double ks_statistic(std::vector<double> &sample1, std::vector<double> &sample2)
{

  std::sort(sample1.begin(), sample1.end());

  std::sort(sample2.begin(), sample2.end());

  double d = 0.0;

  for (size_t i = 0; i < sample1.size(); ++i)
  {
    d = std::max(std::abs(sample1[i] - sample2[i]), d);
  }

  return d;
}

void check(std::vector<double> &theoretical,
           std::vector<double> &samples,
           size_t n_sample)
{
  calculate_correlation(samples, n_sample);

  mean_variance(samples, n_sample);

  double D = ks_statistic(samples, theoretical);

  std::cout << "Kolmogorov-Smirnov : " << D << std::endl;

  double critical_value = 1.36 / sqrt(samples.size() + theoretical.size());
  std::cout << "Critical (5%) : " << critical_value << std::endl;
  if (D < critical_value)
  {
    std::cout << "KS: OK" << std::endl;
  }
  else
  {
    std::cout << "KS: KO" << std::endl;
    ;
  }
}

void fast_uniform_batch()
{
  constexpr size_t n_sample = 10'000;
  std::vector<double> samples(n_sample);
  std::vector<double> theoretical(n_sample);
  std::vector<double> std_sample(n_sample);

  std::mt19937 seed(std::random_device{}());
  std::uniform_real_distribution<> dist(0.0, 1.0);
  MC::KPRNG rng;
  for (size_t i = 0; i < n_sample; ++i)
  {
    samples[i] = rng.double_uniform();
    theoretical[i] = static_cast<double>(i) / n_sample;
    std_sample[i] = dist(seed);
  }

  std::cout << "---REFERENCE---" << std::endl;

  check(theoretical, std_sample, n_sample);

  std::cout << "---MC::RNG---" << std::endl;
  check(theoretical, samples, n_sample);
}

void long_uniform_batch()
{
  constexpr size_t n_sample = 100'000'000;
  std::vector<double> samples(n_sample);
  std::vector<double> theoretical(n_sample);
  std::vector<double> std_sample(n_sample);

  std::mt19937 seed(std::random_device{}());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  MC::KPRNG rng;

  for (size_t i = 0; i < n_sample; ++i)
  {
    samples[i] = rng.double_uniform();
    theoretical[i] = static_cast<double>(i) / n_sample;
    std_sample[i] = dist(seed);
  }

  std::cout << "---REFERENCE---" << std::endl;

  check(theoretical, std_sample, n_sample);

  std::cout << "---MC::RNG---" << std::endl;
  check(theoretical, samples, n_sample);
}

void bench()
{
  std::cout << "BENCH SERIAL:" << std::endl;

  MC::KPRNG rng;
  constexpr size_t n_sample = 100'000'000;

  std::mt19937 seed(std::random_device{}());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  {
    std::vector<double> sample(n_sample);
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < n_sample; ++i)
    {
      sample[i] = dist(seed);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "MC:PRNG duration: " << duration.count() << " milliseconds"
              << std::endl;

    std::cout << "MC::PRNG duration per sample "
              << static_cast<double>(duration.count()) / n_sample
              << " milliseconds" << std::endl;
  }

  {
    std::vector<double> std_sample(n_sample);

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n_sample; ++i)
    {
      std_sample[i] = dist(seed);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "STD duration: " << duration.count() << " milliseconds"
              << std::endl;

    std::cout << "STD duration per sample "
              << static_cast<double>(duration.count()) / n_sample
              << " milliseconds" << std::endl;
  }

  std::cout << "BENCH PARALLEL:" << std::endl;

  {
    Kokkos::Timer timer;
    double time1 = timer.seconds();
    timer.reset();

    auto view = rng.double_uniform(n_sample);

    double time2 = timer.seconds();
    double duration = time2 - time1;
    std::cout << "MC:PRNG duration: " << duration << " milliseconds"
              << std::endl;

    std::cout << "MC::PRNG duration per sample "
              << static_cast<double>(duration) / n_sample << " milliseconds"
              << std::endl;
  }
  std::cout << "BENCH PARRALEL STATIC" << std::endl;
  {

    Kokkos::Timer timer;
    double time1 = timer.seconds();
    timer.reset();

    auto view = rng.double_uniform_view<n_sample>();

    double time2 = timer.seconds();
    double duration = time2 - time1;
    std::cout << "MC:PRNG duration: " << duration << " milliseconds"
              << std::endl;

    std::cout << "MC::PRNG duration per sample "
              << static_cast<double>(duration) / n_sample << " milliseconds"
              << std::endl;
  }
}

int main()
{
  Kokkos::initialize();

  fast_uniform_batch();

  std::cout << "----LONG----" << std::endl;
  long_uniform_batch();

  bench();
  Kokkos::finalize();
}