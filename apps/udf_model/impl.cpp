// #include "mc/particles/data_holder.hpp"
// #include "mc/particles/particle_model.hpp"
// #include "mc/prng/prng.hpp"
// #include <cstdlib>
// #include <dynlib/dyn_module.hpp>
// #include <models/udfmodel_user.hpp>
// #include <models/uptake.hpp>
// #include <models/utils.hpp>
// #include <udf_includes.hpp>

// struct contribs
// {
//   float phi_s;
//   float phi_o;
//   float phi_a;
// };

// namespace Models
// {

//   struct UserImpl
//   {
//     DECLARE_UPTAKE_PROPERTIES(float)
//     float length;
//     float nu1;
//     float nu2;
//     float l_cp;
//     float mu;
//     contribs contrib;
//   };
// } // namespace Models

// namespace
// {

//   constexpr float l_max = 3e-6;
//   constexpr float l_c = 2e-6;
//   constexpr float tau_1 = 1000.;
//   constexpr float tau_2 = 1000.;

//   constexpr float a_11 = 1 / 180.;
//   constexpr float a_21 = 1 / 32.;

//   constexpr float a_12 = 1 / 180.;
//   constexpr float a_22 = 1 / 32.;
//   constexpr float d_m = 0.6e-6;
//   constexpr float minimal_length = 1.4e-6;

//   constexpr float factor = 1000. * 3.14 * d_m * d_m / 4.;
//   constexpr float y_sx = 1. / 2.217737e+00; // y1b

//   float consteval get_phi_s_max(float nu)
//   {
//     return (nu * factor) * y_sx; // kgS/s
//   }

//   constexpr float mt = factor * minimal_length;

//   constexpr float mu_max = 0.8 / 3600.;

//   constexpr float nu_max2 = mu_max * mt;

//   constexpr float nu_max = 8 * 2e-10; // m/s  https://doi.org/10.7554/eLife.67495;

//   constexpr float mu_max_ = nu_max * factor / mt * 3600;

//   constexpr float phi_s_max = get_phi_s_max(nu_max); // 5*mt/7200;
//   constexpr float phi_perm_max = phi_s_max / 40.;    // 5*mt/7200;
//   constexpr float phi_o2_max = 10 * phi_s_max / 180 * 3 * 32;
//   constexpr float NPermease_max = 200;

//   constexpr float y_sxf = y_sx / 3.; // y1
//   // constexpr float y_os = 4.432918e-01;
//   constexpr float y_sa = 0.8;

//   constexpr float y_os = 3; // 3 mol o2 per mol for glucose

//   constexpr float k_pts = 1e-3;
//   constexpr float k_o = 0.0001; // g/L: Anane et. al 2017 (Biochem. Eng. J)

//   constexpr float mmx = 113.1;

//   constexpr float MolarMassG = 180;
//   constexpr float MolarMassO2 = 32;
// } // namespace

// double _mass(Models::UserImpl& pimpl);
// __attribute__((constructor)) static void foo()
// {
//   std::cout << "UDF PIMPL LOADED" << std::endl;
// }

// void _make_udf(Models::User& model)
// {
//   model.pimpl = new Models::UserImpl;
// }

// void _init_udf(Models::UserImpl& pimpl, MC::ParticleDataHolder&, MC::KPRNG _rng)
// {

//   pimpl.nu1 = nu_max * factor / 4;
//   pimpl.nu2 = 0; // pimpl.nu1/5.;

//   pimpl.contrib = {0., 0.0, 0.};
//   auto g = _rng.random_pool.get_state();
//   pimpl.a_pts = g.frand(0.3, 1);
//   pimpl.a_permease = g.frand(0, 1);
//   pimpl.n_permease = g.frand(1, NPermease_max);
//   pimpl.length = l_c; //(float)Kokkos::max((double)minimal_length, g.normal(1e-6, 0.7e-6));
//   pimpl.l_cp =
//       Kokkos::min((float)Kokkos::max((double)minimal_length, g.normal(l_c, l_c / 7.)), l_max);
//   _rng.random_pool.free_state(g);
// }

// float fphi_s(float s)
// {
//   return phi_s_max * s / (1e-3 + s);
// }

// void _update_udf(Models::UserImpl& pimpl,
//                  double d_t,
//                  MC::ParticleDataHolder& p,
//                  const LocalConcentrationView& concentrations)
// {
//   const float s = static_cast<float>(concentrations(0));

//   const float phi_s =
//       Models::Uptake::uptake<float>((float)d_t, pimpl, s, phi_s_max, phi_perm_max, NPermease_max);

//   const float o = concentrations(1);
//   const float phi_o2 = (phi_o2_max)*o / (o + k_o); // gO2/s
//   const float nu_1_star =
//       y_sx * MolarMassG * Kokkos::min(phi_s / MolarMassG, phi_o2 / MolarMassO2 / y_os); // gX/s
//   const float s_1_star = (1 / y_sx * nu_1_star);
//   const float phi_s_residual_1 = Kokkos::max(phi_s - s_1_star, 0.F);
//   const float nu_2_star = y_sxf * phi_s_residual_1;         // gX/s
//   const float nu_eff_1 = Kokkos::min(nu_1_star, pimpl.nu1); // gX/s
//   const float nu_eff_2 = Kokkos::min(nu_2_star, pimpl.nu2); // gX/s
//   const float s_growth = s_1_star + (1 / y_sxf * nu_eff_2);
//   const float s_overflow = phi_s - s_growth;
//   pimpl.mu = (nu_eff_1 + nu_eff_2) / _mass(pimpl);
//   p.status = (pimpl.length > pimpl.l_cp) ? MC::CellStatus::CYTOKINESIS : MC::CellStatus::IDLE;
//   pimpl.contrib.phi_s = -phi_s;
//   pimpl.contrib.phi_o =
//       -1 * ((1. / y_sx / MolarMassG * y_os * MolarMassO2 * nu_eff_1) + 0. * nu_eff_2);
//   pimpl.contrib.phi_a = s_overflow > 0. ? y_sa * (s_overflow) : 0;

//   MSTEP(pimpl.nu1, (nu_1_star - pimpl.nu1) / tau_1);
//   MSTEP(pimpl.nu2, (nu_2_star - pimpl.nu2) / tau_2);
//   MSTEP(pimpl.length, (nu_eff_1 + nu_eff_2) / factor);
// }

// void _contributions_udf(Models::UserImpl& pimpl,
//                         MC::ParticleDataHolder& p,
//                         const ContributionView& contribution)
// {
//   auto access_contribs = contribution.access();

//   access_contribs(0, p.current_container) += p.weight * pimpl.contrib.phi_s;
//   access_contribs(1, p.current_container) += p.weight * pimpl.contrib.phi_o;
//   access_contribs(2, p.current_container) += p.weight * pimpl.contrib.phi_a;
// }

// Models::UserImpl* _division_udf(Models::UserImpl& pimpl, MC::ParticleDataHolder& p, MC::KPRNG rng)
// {
//   const float l = pimpl.length / 2.;

//   auto child_pimpl = new (std::nothrow) Models::UserImpl(pimpl); // NOLINT

//   if (child_pimpl == nullptr)
//   {
//     return nullptr;
//   }
//   pimpl.length = l;
//   child_pimpl->length = l;
//   auto g = rng.random_pool.get_state();
//   child_pimpl->l_cp =
//       Kokkos::min((float)Kokkos::max(minimal_length, (float)g.normal(l_c, l_c / 7.)), l_max);
//   rng.random_pool.free_state(g);
//   return child_pimpl;
// }

// double _mass(Models::UserImpl& pimpl)
// {
//   return pimpl.length * factor;
// }

// void _fill_properties(Models::UserImpl& pimpl, SubViewtype full)
// {
//   full(0) = _mass(pimpl);
//   full(1) = pimpl.length * 1e6;
//   full(2) = pimpl.nu1;
//   full(3) = pimpl.nu2;
//   full(4) = pimpl.mu;
//   full(5) = pimpl.a_pts;
//   full(6) = pimpl.a_permease;
// }

// std::vector<std::string> _names()
// {
//   return {"mass", "length", "nu1", "nu2", "mu", "a_pts", "a_permease"};
// }

// size_t _number()
// {
//   return _names().size();
// }

// //

// DECLARE_DELETER(Models::UserImpl)

// EXPORT_MODULE(module,
//               &_init_udf,
//               &_make_udf,
//               &_delete_udf,
//               &_update_udf,
//               &_contributions_udf,
//               &_division_udf,
//               &_mass,
//               &_fill_properties,
//               &_names,
//               &_number);
