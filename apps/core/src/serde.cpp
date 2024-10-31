#include "mc/unit.hpp"
#include "simulation/scalar_initializer.hpp"
#include "simulation/simulation.hpp"
#include <cstddef>
#include <ios>
#include <memory>
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
#  include <cereal/types/vector.hpp> //MC::List use vector internally
#  include <fstream>
#  include <optional>
#  include <serde.hpp>
#  include <string_view>

static void write_to_file(const std::ostringstream &oss, std::string_view filename)
{

  std::ofstream file(filename.data());
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
}

#  define PLAINTEXT_ARCHIVE cereal::XMLOutputArchive
#  define BINARY_ARCHIVE cereal::BinaryArchive

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

  void save_simulation(const Core::CaseData &case_data)
  {

    std::stringstream serde_name;
    serde_name << case_data.params.user_params.results_file_name << "_serde_" << case_data.exec_info.current_rank;

    std::ostringstream buf(std::ios::binary);
    {

      Archive_t ar(buf);
      ar(cereal::make_nvp("version", ExecInfo::get_version()), case_data.exec_info);
      auto dim = case_data.simulation->getDimensions();
      auto cliq = case_data.simulation->getCliqData();
      auto cgas = case_data.simulation->getCgasData();

      std::optional<std::vector<double>> cgas_a =
          cgas.has_value() ? std::make_optional(std::vector<double>(cgas->begin(), cgas->end())) : std::nullopt;

      ar(dim, std::vector<double>(cliq.begin(), cliq.end()), cgas_a);

      ar(cereal::defer(case_data.simulation->mc_unit));
      ar.serializeDeferments();
    }

    write_to_file(buf, serde_name.str());
  }

  bool load_simulation(Core::CaseData& res,std::string_view ser_filename)
  {

    std::stringstream buffer;
    read_file(buffer, ser_filename);
    iArchive_t ar(buffer);

    std::string version;
    ar(version, res.exec_info);

    std::cout << res.exec_info.run_id << std::endl;

    Simulation::Dimensions dims;
    std::vector<double> cliq;
    std::optional<std::vector<double>> cgas;

    ar(dims, cliq, cgas);

    Simulation::ScalarInitializer scalar_init;
    scalar_init.n_species = dims.n_species;
    scalar_init.type = Simulation::ScalarInitialiserType::File;
    scalar_init.liquid_buffer = cliq;
    scalar_init.gas_buffer = cgas;



    if(cgas)
    {
      std::cout<<cgas->size()<<"\r\n";
    }

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;
    ar(mc_unit);
    res.simulation =std::make_unique<Simulation::SimulationUnit>(std::move(mc_unit),scalar_init);

    return true;
  }

} // namespace SerDe

#endif // USE_CEAREAL