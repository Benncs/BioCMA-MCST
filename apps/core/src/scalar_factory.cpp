
#include <algorithm>
#include <cma_read/light_2d_view.hpp>
#include <core/scalar_factory.hpp>
#include <filesystem>
#include <optional>
#include <simulation/scalar_initializer.hpp>
#include <span>
#include <stdexcept>
#include <utility>
#include <variant>
#ifdef USE_HIGHFIVE

#ifndef NDEBUG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif 
#  include <Eigen/Dense>
#ifndef NDEBUG
#pragma GCC diagnostic pop
#endif 
#  include <highfive/H5DataSet.hpp>
#  include <highfive/H5DataSpace.hpp>
#  include <highfive/H5Easy.hpp>
#  include <highfive/H5File.hpp>
#  include <highfive/H5PropertyList.hpp>
#endif

namespace Core::ScalarFactory
{

  Simulation::ScalarInitializer scalar_factory(bool /*f_init_gas_flow*/,
                                               std::span<double> gas_volume,
                                               std::span<double> liquid_volume,
                                               ScalarVariant arg_liq)
  {
    if (gas_volume.size() != liquid_volume.size())
    {
      throw std::invalid_argument("Error size of volume need to be the same");
    }

    auto ret = std::visit(Visitor(), std::move(arg_liq));

    ret.volumesgas = gas_volume;
    ret.volumesliq = liquid_volume;

    // FIXME
    //  if (ret.gas_flow != f_init_gas_flow)
    //  {
    //    throw std::invalid_argument("Gas provided but no functor, or inverse");
    //  }

    if (!sanitize(ret))
    {
      throw std::invalid_argument("Error while initializing");
    }

    return ret;
  }

  Simulation::ScalarInitializer Visitor::operator()(Uniform args)
  {
    auto res = Simulation::ScalarInitializer();
    res.type = Simulation::ScalarInitialiserType::Uniform;
    const size_t n_species = args.liquid_concentrations.size();
    res.n_species = n_species;
    res.gas_flow = false;

    auto wrap_functor = [](auto&& c)
    {
      // return [n_species, local_concentrations = c](size_t i, CmaRead::L2DView<double>& view)
      // {
      //   for (size_t i_species = 0; i_species < n_species; ++i_species)
      //   {
      //     view(i_species, i) = local_concentrations[i_species];
      //   };
      // };
      return [local_concentrations = c](size_t i_row, size_t i_col)
      {
        (void)i_col;
        return local_concentrations[i_row];
      };
    };

    res.liquid_f_init = wrap_functor(args.liquid_concentrations);

    if (args.gas_concentration.has_value())
    {
      res.gas_f_init = wrap_functor(*args.gas_concentration);
      res.gas_flow = true;
    }

    return res;
  }

  // Simulation::ScalarInitializer Visitor::operator()(CustomFunctor func)
  // {
  //   auto res = Simulation::ScalarInitializer();
  //   res.type = Simulation::ScalarInitialiserType::CustomFunctor;
  //   res.n_species = func.n_species;

  //   res.gas_flow = false;

  //   res.liquid_f_init = func.liquid;

  //   if (func.gas.has_value())
  //   {
  //     res.gas_f_init = *func.gas;
  //     res.gas_flow = true;
  //   }

  //   return res;
  // }

  Simulation::ScalarInitializer Visitor::operator()(Local args)
  {
    auto res = Simulation::ScalarInitializer();
    res.type = Simulation::ScalarInitialiserType::Local;
    const auto indices = args.liquid_indices;
    const auto concentrations = args.liquid_concentrations;
    const size_t n_species = concentrations.size();
    res.n_species = n_species;
    res.gas_flow = false;

    // const auto wrap_functor = [n_species](auto&& i, auto&& c)
    // {
    //   return [n_species,
    //           concentrations = std::forward<decltype(c)>(c),
    //           _indices = std::forward<decltype(i)>(i)](size_t i, CmaRead::L2DView<double>& view)
    //   {
    //     if (std::ranges::find(_indices, i) != _indices.end())
    //     {
    //       for (size_t i_species = 0; i_species < n_species; ++i_species)
    //       {
    //         view(i_species, i) = concentrations[i_species];
    //       }
    //     }
    //   };
    // };

    const auto wrap_functor = [](auto&& i, auto&& c)
    {
      return [concentrations = std::forward<decltype(c)>(c),
              _indices = std::forward<decltype(i)>(i)](std::size_t i_row, std::size_t i_col)
      {
        if (std::ranges::find(_indices, i_col) != _indices.end())
        {
          return concentrations[i_row];
        }
        else
        {
          return 0.;
        }
      };
    };

    res.liquid_f_init = wrap_functor(indices, concentrations);

    if (args.gas_concentration.has_value() && args.gas_indices.has_value())
    {
      res.gas_f_init = wrap_functor(*args.gas_indices, *args.gas_concentration);
      res.gas_flow = true;
    }

    return res;
  }

  Simulation::ScalarInitializer Visitor::operator()(FullCase data)
  {
    auto res = Simulation::ScalarInitializer();
    res.type = Simulation::ScalarInitialiserType::FullCase;
    res.gas_flow = data.raw_gas.has_value();
    res.n_species = data.n_species;
    if (res.gas_flow)
    {
      if (data.raw_gas->size() != data.raw_liquid.size())
      {
        throw std::invalid_argument("Liquid and Gas data should have the same size");
      }
      res.gas_buffer = std::move(*data.raw_gas);
    }

    res.liquid_buffer = std::move(data.raw_liquid);

    return res;
  }

  Simulation::ScalarInitializer Visitor::operator()(File filepath)
  {

#ifdef USE_HIGHFIVE
    auto res = Simulation::ScalarInitializer();
    res.type = Simulation::ScalarInitialiserType::File;
    if (!std::filesystem::is_regular_file(filepath.path))
    {
      throw std::invalid_argument("Unable to open provided concentration initaliser file");
    }

    HighFive::File file(filepath.path.data(), HighFive::File::ReadOnly);

    auto liquid_dataset = file.getDataSet("initial_liquid");

    HighFive::DataSet gas_dataset;
    bool gas = false;
    if (file.exist("initial_gas"))
    {
      gas_dataset = file.getDataSet("initial_gas");
      gas = true;
    }

    auto set_buffer = [&filepath, &res](auto& dataset)
    {
      auto dims = dataset.getDimensions();

      auto n_elements = dims[0] * dims[1];
      if (n_elements >= filepath.n_compartment)
      {
        auto nd_array = std::vector<double>(n_elements);
        dataset.template read_raw<double>(nd_array.data());
        res.n_species = n_elements / filepath.n_compartment;
        return nd_array;
      }

      throw std::invalid_argument(__FILE__ "Invalid file, size don't match");
    };

    res.liquid_buffer = set_buffer(liquid_dataset);
    if (gas)
    {
      res.gas_buffer = set_buffer(gas_dataset);
      res.gas_flow = true;
    }

    // First read the dimensions.
    // auto dims = liquid_dataset.getDimensions();

    // auto n_elements = dims[0] * dims[1];
    // if (n_elements >= arg.n_compartment)
    // {
    //   auto nd_array = std::vector<double>(n_elements);
    //   liquid_dataset.read_raw<double>(nd_array.data());
    //   res.liquid_buffer = std::move(nd_array);
    //   res.n_species = n_elements / arg.n_compartment;
    // }

    return res;
#else
    throw std::invalid_argument("Not implemented yet");
#endif
  }

  Simulation::ScalarInitializer Visitor::operator()(CustomScript /*path*/)
  {
    throw std::invalid_argument("Not implemented yet");
  }

  bool sanitize(const Simulation::ScalarInitializer& res)
  {

    bool flag = false;

    auto test_functor = [](auto&& _res)
    {
      bool _flag = _res.liquid_f_init.has_value() && !_res.liquid_buffer.has_value();

      if (_res.gas_flow)
      {
        _flag = _flag && _res.gas_f_init.has_value();
      }
      return _flag;
    };
    switch (res.type)
    {

    case Simulation::ScalarInitialiserType::Uniform: // NOLINT
    {
      flag = test_functor(res);
      break;
    }

    case Simulation::ScalarInitialiserType::Local:
    {
      flag = test_functor(res);
      break;
    }
    case Simulation::ScalarInitialiserType::File:
    {
      // File doesn't need functor but buffer
      // First check initialiser has buffer and functor and not set
      flag = res.liquid_buffer.has_value() &&
             (!res.gas_f_init.has_value() && !res.liquid_f_init.has_value());
      if (flag)
      {
        flag = res.liquid_buffer->size() != 0; // If buffer check that is not empty
      }
      break;
    }
    case Simulation::ScalarInitialiserType::CustomScript:
    {
      flag = false;
      break;
    }
    case Simulation::ScalarInitialiserType::FullCase:
    {
      flag = true; // TODO
      break;
    }
    }

    return flag;
  }

} // namespace Core::ScalarFactory
