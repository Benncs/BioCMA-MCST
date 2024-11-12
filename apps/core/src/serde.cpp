#include "mc/unit.hpp"
#include "simulation/scalar_initializer.hpp"
#include "simulation/simulation.hpp"
#include <cstddef>
#include <exception>
#include <ios>
#include <memory>
#include <stdexcept>
#include <vector>
#ifdef USE_CEAREAL
#  include "common/execinfo.hpp"
#  include <cereal/archives/binary.hpp>
#  include <cereal/archives/xml.hpp>
#  include <cereal/types/array.hpp> //MC::events use vector internally
#  include <cereal/types/memory.hpp>
#  include <cereal/types/optional.hpp>
#  include <cereal/types/string.hpp> //MC::List use vector internally
#  include <cereal/types/tuple.hpp>
#  include <cereal/types/variant.hpp> //MC::Unit use vector internally
#  include <cereal/types/vector.hpp>  //MC::List use vector internally
#  include <fstream>
#  include <optional>
#  include <serde.hpp>
#  include <string_view>

static void write_to_file(const std::ostringstream &oss, std::string_view filename)
{

  std::ofstream file(filename.data(), std::ios::binary);
  if (!file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  file << oss.str();
  file.close();
}

static void read_file(std::stringstream &buffer, std::string_view filename)
{
  std::fstream file;
  file.open(filename.data(), std::ios::binary | std::ios::in);
  if (file.is_open())
  {
    buffer << file.rdbuf();
    file.close();
  }
  else
  {
    throw std::runtime_error("cannot read file");
  }
}

using Archive_t = cereal::XMLOutputArchive;
using iArchive_t = cereal::XMLInputArchive;

// void read_archive(cereal::XMLInputArchive& ar,std::string_view filename)
// {
//   std::fstream file;
//   file.open(filename.data(), std::ios::binary | std::ios::in);
//   std::stringstream buffer;
//   if (file.is_open())
//   {
//     buffer << file.rdbuf();
//     file.close();
//   }
//   ar(buffer);
//   file.close();
// };

// void read_archive(cereal::BinaryInputArchive& ar,std::string_view filename)
// {

// };

namespace SerDe
{
  static std::vector<double> spl = {20e-3};
  static std::vector<double> spg = {};
  void save_simulation(const Core::CaseData &case_data)
  {

    std::stringstream serde_name;
    serde_name << case_data.params.results_file_name << "_serde_" << case_data.exec_info.current_rank
               << ".raw";

    std::ostringstream buf(std::ios::binary);
    {

      Archive_t ar(buf);
      ar(cereal::make_nvp("version", ExecInfo::get_version()), case_data.exec_info);
      auto dim = case_data.simulation->getDimensions();
      auto cliq = case_data.simulation->getCliqData();
      auto cgas = case_data.simulation->getCgasData();

      std::optional<std::vector<double>> cgas_a =
          cgas.has_value() ? std::make_optional(std::vector<double>(cgas->begin(), cgas->end())) : std::nullopt;

      ar(dim, std::vector<double>(cliq.begin(), cliq.end()), cgas_a, case_data.simulation->get_end_time_mut());
      std::cerr<<"DOMAIN"<<case_data.simulation->mc_unit->domain.getNumberCompartments()<<std::endl;
      ar(case_data.simulation->mc_unit);
    }

    write_to_file(buf, serde_name.str());
  }

  bool load_simulation(Core::GlobalInitialiser &gi, Core::CaseData &case_data, std::string_view ser_filename)
  {

    std::stringstream buffer;

    read_file(buffer, ser_filename);

    iArchive_t ar(buffer);

    std::string version;

    ExecInfo serde_exec{};
    ar(version, serde_exec);

    case_data.exec_info.run_id = serde_exec.run_id;
    std::cout << "SIMULATION: " << case_data.exec_info.run_id << " LOADED" << std::endl;

    Simulation::Dimensions dims;
    std::vector<double> read_c_liq;
    std::optional<std::vector<double>> read_c_gas;
    double start_time{};
    ar(dims, read_c_liq, read_c_gas, start_time);
    
    // const auto *state = case_data.transitioner->getState();

    Simulation::ScalarInitializer scalar_init;
    scalar_init.n_species = dims.n_species;
    scalar_init.type = Simulation::ScalarInitialiserType::File;
    scalar_init.liquid_buffer = read_c_liq;
    scalar_init.gas_buffer = read_c_gas;
    scalar_init.gas_flow = read_c_gas.has_value();

    // FIXME
    scalar_init.volumesliq = spl;
    scalar_init.volumesgas = spg;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;
    ar(mc_unit);
    std::cerr<<"DOMAIN"<<mc_unit->domain.getNumberCompartments()<<std::endl;
    auto simulation = gi.init_simulation(std::move(mc_unit), scalar_init);

    if (!simulation.has_value())
    {
      return false;
    }

    case_data.simulation = std::move(*simulation);

    case_data.simulation->get_start_time_mut() = start_time;
    return true;
  }

} // namespace SerDe

#endif // USE_CEAREAL