#include "mc/particles/particle_model.hpp"
#include "mc/prng/prng.hpp"
#include <cstdlib>
#include <dynlib/dyn_module.hpp>
#include <models/model_user.hpp>
#include <udf_includes.hpp>
namespace Models
{
  struct UserImpl
  {
    float a;
    
  };
} // namespace Models

// __attribute__((constructor)) static void foo()
// {
//   std::cout << "IMPL LOADED" << std::endl;
// }

void _init_udf(Models::UserImpl& pimpl, MC::ParticleDataHolder&)
{
}

void _make_udf(Models::User& model)
{
  model.pimpl = new Models::UserImpl;
}

void _update_udf(Models::UserImpl& pimpl, double d_t, MC::ParticleDataHolder& p,const LocalConcentrationView& concentrations)
{
  pimpl.a = concentrations[0];
}

void _contributions_udf(Models::UserImpl& pimpl, MC::ParticleDataHolder& p,
                                          ContributionView& contribution)
{
  auto access_contribs = contribution.access();
  
  double ra = 0.1e-15*pimpl.a*pimpl.a;

  access_contribs(0,p.current_container)+= (-ra*p.weight);
  access_contribs(1,p.current_container)+= (2*ra*p.weight);
  
}

Models::UserImpl* _division_udf(Models::UserImpl& pimpl, MC::ParticleDataHolder& p,
                                          MC::KPRNG )
{
  return new Models::UserImpl(pimpl);
}

DECLARE_DELETER(Models::UserImpl)

EXPORT_MODULE(module, &_init_udf, &_make_udf, &_delete_udf, &_update_udf,&_contributions_udf,&_division_udf);
