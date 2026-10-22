#include <Kokkos_Core.hpp>
#include <cassert>
#include <common/common.hpp>
#include <common/traits.hpp>
#include <cstddef>
#include <hydro/impl_mass_transfer.hpp>
#include <memory>
#include <mixture/species_descriptor.hpp>
#include <scalar.hpp>
#include <simulation/mass_transfer.hpp>
#include <stdexcept>
#include <vector>

constexpr double TEST_TOLERANCE = 1e-12;

constexpr double o2_molar_weight = 32.;
constexpr double o2_henry = 3e-2;
constexpr double glucose_molar_weight = 180.;
constexpr double kla = 0.5;
constexpr double volume = 2.;

std::shared_ptr<Simulation::impl::ScalarSimulation>
mock_scalar(std::size_t n_species, std::size_t n_compartments, double volume)
{
  std::vector<double> volumes(n_compartments, volume);
  return std::make_shared<Simulation::impl::ScalarSimulation>(
      n_compartments, n_species, volumes, volumes);
}

void
test_mtr_fixed_kla()
{
  constexpr std::size_t n_species = 1;
  constexpr std::size_t n_compartments = 2;
  constexpr double liquid_0 = 1e-3;
  constexpr double liquid_1 = 2e-3;
  constexpr double gas_0 = 0.3;
  constexpr double gas_1 = 0.2;

  auto table = Mixture::SpecieTable{ Mixture::new_specie(
      "o2", o2_molar_weight, o2_henry) };

  auto scalar = mock_scalar(n_species, n_compartments, volume);
  auto liquid = scalar->concentration_data(Phase::Liquid);
  auto gas = scalar->concentration_data(Phase::Gas);
  liquid[0] = liquid_0;
  liquid[1] = liquid_1;
  gas[0] = gas_0;
  gas[1] = gas_1;

  Simulation::MassTransfer::Type::FixedKla fixed_kla;
  fixed_kla.value = { kla };

  auto model = Simulation::MassTransfer::MassTransferModel(
      table, fixed_kla, scalar->mtr_input());

  model.gas_liquid_mass_transfer();

  const auto& mtr = model.proxy()->mtr;
  for (std::size_t i = 0; i < n_compartments; ++i)
  {
    auto expected = kla * (o2_henry * gas[i] - liquid[i])
                    * scalar->volume_at(Phase::Liquid, i);
    assert(almost_equal(mtr(0, EIGEN_INDEX(i)), expected, TEST_TOLERANCE));
  }
}

void
test_mtr_no_henry()
{
  constexpr std::size_t n_species = 1;
  constexpr std::size_t n_compartments = 1;
  constexpr double liquid_0 = 5.;
  constexpr double gas_0 = 1.;

  auto table = Mixture::SpecieTable{ Mixture::new_specie(
      "glucose", glucose_molar_weight, std::nullopt) };

  auto scalar = mock_scalar(n_species, n_compartments, volume);
  scalar->concentration_data(Phase::Liquid)[0] = liquid_0;
  scalar->concentration_data(Phase::Gas)[0] = gas_0;

  Simulation::MassTransfer::Type::FixedKla fixed_kla;
  fixed_kla.value = { kla };

  auto model = Simulation::MassTransfer::MassTransferModel(
      table, fixed_kla, scalar->mtr_input());

  model.gas_liquid_mass_transfer();

  assert(model.proxy()->mtr(0, 0) == 0.);
}

void
test_mtr_kla_size()
{
  constexpr std::size_t n_species = 1;
  constexpr std::size_t n_compartments = 1;

  auto table = Mixture::SpecieTable{ Mixture::new_specie(
      "o2", o2_molar_weight, o2_henry) };

  auto scalar = mock_scalar(n_species, n_compartments, volume);

  Simulation::MassTransfer::Type::FixedKla fixed_kla;
  fixed_kla.value = { kla, kla };

  bool thrown = false;
  try
  {
    auto model = Simulation::MassTransfer::MassTransferModel(
        table, fixed_kla, scalar->mtr_input());
  }
  catch (const std::invalid_argument&)
  {
    thrown = true;
  }
  assert(thrown);
}

int
main()
{
  Kokkos::initialize();
  test_mtr_fixed_kla();
  test_mtr_no_henry();
  test_mtr_kla_size();
  Kokkos::finalize();
  return 0;
}
