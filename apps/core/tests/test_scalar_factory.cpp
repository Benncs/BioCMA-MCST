#include <cassert>
#include <filesystem>
#include <iostream>
#include <core/scalar_factory.hpp>
#include <stdexcept>
#include <string_view>
#include <vector> 

#define WRAP_EXCEP(__f__)                                                      \
  try                                                                          \
  {                                                                            \
    __f__;                                                                     \
  }                                                                            \
  catch (std::exception & e)                                                   \
  {                                                                            \
    std::cout << __LINE__ << std::endl;                                        \
    std::cout << e.what() << std::endl;                                        \
    assert(false);                                                             \
  }

static constexpr size_t n_compartment = 10;
static constexpr size_t n_species = 3;
#define GET_VOLUME                                                             \
  static std::vector<double> gas_volume = {0., 0., 0., 0.};                    \
  static std::vector<double> liquid_volume = {1., 1., 1., 1.};

#define GET_CONCENTRATION std::vector<double> concentrations = {1., 5., 6.};

#define GET_CONCENTRATION_GAS                                                  \
  std::vector<double> concentrations_gas = {0.1, 0.2, 3.};

static std::vector<double> get_raw_concentration_data()
{
  return std::vector<double>(n_compartment * n_species);
}

void mock_init(CmaRead::L2DView<double> c, auto functor)
{
  assert(functor.has_value());
  for (size_t i = 0; i < c.getNCol(); ++i)
  {
    (*functor)(i, c);
  }
}

#define DO_INIT_LIQ                                                            \
  auto rd = get_raw_concentration_data();                                      \
  auto cliq = CmaRead::L2DView<double>(rd, n_species, n_compartment);          \
  mock_init(cliq, scalar_init.liquid_f_init);

#define DO_INIT_GAS                                                            \
  auto rdg = get_raw_concentration_data();                                     \
  auto cgas = CmaRead::L2DView<double>(rdg, n_species, n_compartment);         \
  mock_init(cgas, scalar_init.gas_f_init);

void test_uniform_liq()
{
  GET_VOLUME
  GET_CONCENTRATION

  bool f_init_gas_flow = false;

  Core::ScalarFactory::Uniform arg = {concentrations, std::nullopt};

  auto scalar_init = Core::ScalarFactory::scalar_factory(
      f_init_gas_flow, gas_volume, liquid_volume, arg);

  assert(scalar_init.n_species == n_species);

  DO_INIT_LIQ

  for (size_t i_compartment = 0; i_compartment < n_compartment; ++i_compartment)
  {
    for (size_t isp = 0; isp < n_species; ++isp)
    {
      assert(cliq(isp, i_compartment) == concentrations[isp]);
    }
  }
}

void test_uniform_liq_gas()
{
  GET_VOLUME
  GET_CONCENTRATION
  GET_CONCENTRATION_GAS
  bool f_init_gas_flow = true;
  Core::ScalarFactory::Uniform arg = {concentrations, concentrations_gas};

  auto scalar_init = Core::ScalarFactory::scalar_factory(
      f_init_gas_flow, gas_volume, liquid_volume, std::move(arg));

  assert(scalar_init.n_species == n_species);

  DO_INIT_LIQ
  DO_INIT_GAS

  for (size_t i_compartment = 0; i_compartment < n_compartment; ++i_compartment)
  {
    for (size_t isp = 0; isp < n_species; ++isp)
    {
      assert(cliq(isp, i_compartment) == concentrations[isp]);
      assert(cgas(isp, i_compartment) == concentrations_gas[isp]);
    }
  }
}

void test_local_gas_liq()
{

  auto test_f = [](size_t __ic, auto &&i, auto &&c, auto &v)
  {
    if (__ic == i[0] || __ic == i[1] || __ic == i[2])
    {
      for (size_t isp = 0; isp < n_species; ++isp)
      {
        assert(v(isp, __ic) == c[isp]);
      }
    }
    else
    {
      for (size_t isp = 0; isp < n_species; ++isp)
      {

        assert(v(isp, __ic) == 0);
      }
    }
  };

  GET_VOLUME GET_CONCENTRATION GET_CONCENTRATION_GAS std::vector<size_t>
      liq_indices = {0, 5, 6};

  std::vector<size_t> gas_indices = {6, 5, 7};

  bool f_init_gas_flow = true;

  Core::ScalarFactory::Local arg = {
      concentrations, liq_indices, concentrations_gas, gas_indices};

  auto scalar_init =
      scalar_factory(f_init_gas_flow, gas_volume, liquid_volume, arg);

  DO_INIT_LIQ
  DO_INIT_GAS

  for (size_t i_compartment = 0; i_compartment < n_compartment; ++i_compartment)
  {
    test_f(i_compartment, liq_indices, concentrations, cliq);
    test_f(i_compartment, gas_indices, concentrations_gas, cgas);
  }
}

void test_local_liq()
{
  GET_VOLUME
  GET_CONCENTRATION

  std::vector<size_t> indices = {0, 5, 6};

  bool f_init_gas_flow = false;

  Core::ScalarFactory::Local arg = {concentrations, indices};

  auto scalar_init =
      scalar_factory(f_init_gas_flow, gas_volume, liquid_volume, arg);

  DO_INIT_LIQ

  for (size_t i_compartment = 0; i_compartment < n_compartment; ++i_compartment)
  {

    if (i_compartment == indices[0] || i_compartment == indices[1] ||
        i_compartment == indices[2])
    {
      for (size_t isp = 0; isp < n_species; ++isp)
      {
        assert(cliq(isp, i_compartment) == concentrations[isp]);
      }
    }
    else
    {
      for (size_t isp = 0; isp < n_species; ++isp)
      {
        std::cerr << cliq(isp, i_compartment) << std::endl;
        ;
        assert(cliq(isp, i_compartment) == 0);
      }
    }
  }
}

void test_read(std::string_view path)
{
  GET_VOLUME
  GET_CONCENTRATION

  auto args = Core::ScalarFactory::File(500,  path);
  Core::ScalarFactory::scalar_factory(false, gas_volume, liquid_volume, args);
}

void test_wrong_size()
{
  GET_CONCENTRATION
  GET_VOLUME
  gas_volume.emplace_back(1); // So that gas and liq don't have the same size
  Core::ScalarFactory::Uniform arg = {concentrations};

  // Should throw exception
  try
  {
    auto scalar_init = scalar_factory(false, gas_volume, liquid_volume, arg);
  }
  catch (std::invalid_argument &e)
  {
    return;
  }
  assert(false);
}

int main(int argc, char **argv)
{

  if (argc == 2)
  {
    std::string_view path = argv[1];
    if (std::filesystem::is_regular_file(path))
    {
      WRAP_EXCEP(test_read(path));
    }
  }

  WRAP_EXCEP(test_wrong_size());
  WRAP_EXCEP(test_uniform_liq());
  WRAP_EXCEP(test_uniform_liq_gas());
  WRAP_EXCEP(test_local_liq());
  WRAP_EXCEP(test_local_gas_liq());

  return 0;
}