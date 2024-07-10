#ifndef __MC_PRNG_HPP__
#define __MC_PRNG_HPP__

#include <cstdint>
#include <mc/prng/distribution.hpp>
#include <random>

inline unsigned tau_step(unsigned &z, int S1, int S2, int S3, unsigned M)
{
  unsigned b = (((z << S1) ^ z) >> S2);
  return z = (((z & M) << S3) ^ b);
}

inline unsigned LCGStep(unsigned &z, unsigned A, unsigned C)
{
  return z = (A * z + C);
}
namespace MC
{

  static constexpr size_t MC_RAND_DEFAULT_SEED = 10;
  class PRNG
  {
  public:
    explicit PRNG()
    {
      gen = std::mt19937(std::random_device{}());
    }

    explicit PRNG(uint64_t seed)
    {
      gen = std::mt19937(seed);
    }

    double uniform_double_rand(double min, double max)
    {
      double rn = 0;
      std::uniform_real_distribution<double> double_distribution(min, max);
      rn = double_distribution(gen);
      return rn;
    }

    inline double double_unfiform()
    {
      return step();
      // return _uniform_double(gen);
    }

    inline float step()
    {
      // Combined period is lcm(p1,p2,p3,p4)~ 2^121
      // https://indico.cern.ch/event/93877/contributions/2118070/attachments/1104200/1575343/acat3_revised_final.pdf
      return 2.3283064365387e-10 *
             (                                            // Periods
                 tau_step(z1, 13, 19, 12, 4294967294UL) ^ // p1=2^31-1
                 tau_step(z2, 2, 25, 4, 4294967288UL) ^   // p2=2^30-1
                 tau_step(z3, 3, 11, 17, 4294967280UL) ^  // p3=2^28-1
                 LCGStep(z4, 1664525, 1013904223UL)       // p4=2^32
             );
    }

    auto &rng()
    {
      return gen;
    };

    static inline auto get_rng(uint64_t seed = MC_RAND_DEFAULT_SEED)
    {
      std::mt19937 _gen(seed);
      return _gen;
    };

  private:
    unsigned z2 = 1;
    unsigned z3 = 2;
    unsigned z4 = 4;
    std::mt19937 gen;
    unsigned z1 = 0;

    std::uniform_real_distribution<double> _uniform_double =
        std::uniform_real_distribution<double>(0., 1.);
  };

} // namespace MC

#endif //__MC_PRNG_HPP__