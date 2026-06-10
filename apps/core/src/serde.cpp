

#ifdef USE_CEAREAL
#  include <cereal/archives/binary.hpp>
#  include <cereal/archives/xml.hpp>
#  include <cereal/types/array.hpp> //MC::events use array internally
#  include <cereal/types/memory.hpp>
#  include <cereal/types/optional.hpp>
#  include <cereal/types/string.hpp> //MC::List use vector internally
#  include <cereal/types/tuple.hpp>
#  include <cereal/types/variant.hpp> //MC::Unit use variant internally
#  include <cereal/types/vector.hpp>  //MC::List use vector internally
#  include <common/execinfo.hpp>
#  include <core/global_initaliser.hpp>
#  include <cstdint>
#  include <fstream>
#  include <ios>
#  include <mc/unit.hpp>
#  include <memory>
#  include <optional>
#  include <serde.hpp>
#  include <simulation/scalar_initializer.hpp>
#  include <simulation/simulation.hpp>
#  include <stdexcept>
#  include <string_view>
#  include <vector>
static void
write_to_file(const std::ostringstream& oss, std::string_view filename)
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

static void
read_file(std::stringstream& buffer, std::string_view filename)
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

// using Archive_t = cereal::XMLOutputArchive;
// using iArchive_t = cereal::XMLInputArchive;
using Archive_t = cereal::BinaryOutputArchive;
using iArchive_t = cereal::BinaryInputArchive;

namespace SerDe
{

  void
  save_simulation(const Core::CaseData& case_data)
  {

    std::stringstream serde_name;
    serde_name << case_data.params.results_file_name << "_serde_"
               << case_data.exec_info.current_rank << ".raw";

    std::ostringstream buf(std::ios::binary);
    {
      auto accessor = case_data.simulation->getter();
      Archive_t ar(buf);
      ar(cereal::make_nvp("version", ExecInfo::get_version()),
         case_data.exec_info);
      auto dim = accessor.getDimensions();
      auto cliq = accessor.getCliqData();

      if (!case_data.simulation->checkScalar())
      {
        throw std::runtime_error(
            "Simulation with negative value won´t be able to be loaded");
      }

      auto cgas = accessor.getCgasData();

      std::optional<std::vector<double>> cgas_a
          = cgas.has_value() ? std::make_optional(std::vector<double>(
                                   cgas->begin(), cgas->end()))
                             : std::nullopt;

      ar(case_data.params.number_particle,
         dim,
         std::vector<double>(cliq.begin(), cliq.end()),
         cgas_a,
         accessor.absolute_time());

      ar(accessor.mc_unit());
    }
    write_to_file(buf, serde_name.str());
  }

  std::optional<Simulation::ScalarInitializer>
  build_scalar_init(Core::GlobalInitialiser& gi,
                    const Simulation::Dimensions& dims,
                    std::vector<double>&& c_liq_buffer,
                    std::optional<std::vector<double>>&& c_gas_buffer)
  {
    auto sc = gi.init_scalar();

    if (!sc.has_value())
    {

      return std::nullopt;
    }

    if (c_liq_buffer.size() % sc->volumesliq.size() != 0)
    {
      throw std::invalid_argument(
          "Liquid Concentration buffer and CM volume have incompatible sizes ");
    }

    if (c_gas_buffer)
    {
      if (c_gas_buffer->size() % sc->volumesgas.size() != 0)
      {
        throw std::invalid_argument(
            "Gas Concentration buffer and CM volume have incompatible sizes ");
      }
    }

    // FIXME
    // Overwrite value set by initialiser
    sc->liquid_f_init = std::nullopt;
    sc->gas_f_init = std::nullopt;
    sc->gas_buffer = std::move(c_gas_buffer);
    sc->liquid_buffer = std::move(c_liq_buffer);
    sc->n_species = dims.n_species;
    sc->type = Simulation::ScalarInitialiserType::Serde;

    return sc;
  }

  bool
  load_simulation(Core::GlobalInitialiser& gi,
                  Core::CaseData& case_data,
                  std::string_view ser_filename)
  {

    std::stringstream buffer;
    read_file(buffer, ser_filename);

    iArchive_t ar(buffer);

    std::string version;
    ExecInfo serde_exec{};
    ar(version, serde_exec);

    case_data.exec_info.run_id = serde_exec.run_id;

    uint64_t np = 0;
    Simulation::Dimensions dims;
    std::vector<double> read_c_liq;
    std::optional<std::vector<double>> read_c_gas;
    double start_time{};
    ar(np, dims, read_c_liq, read_c_gas, start_time);

    auto sc = build_scalar_init(
        gi, dims, std::move(read_c_liq), std::move(read_c_gas));

    if (!sc.has_value())
    {
      std::cout << "Scalar loaded failed" << std::endl;
      return false;
    }

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;
    ar(mc_unit);
    assert(mc_unit != nullptr);

#  warning message("MTR model is not loaded")
    auto simulation = gi.init_simulation(std::move(mc_unit), std::move(*sc));

    if (!simulation.has_value())
    {
      std::cout << "SIMULATION loaded failed" << std::endl;
      return false;
    }
    // TODO the should be useless because of checks inside the simulation
    // constructor
    if (!(*simulation)->checkScalar())
    {
      std::cout << "SIMULATION loaded failed" << std::endl;
      return false;
    }

    case_data.simulation = std::move(*simulation);
    std::vector<double> kla(
        case_data.simulation->getter().getDimensions().n_species);
    if (kla.size() > 1)
    {
      kla[1] = 0.2; // 700 h-1
    }

    auto auto_mtr_type = Simulation::MassTransfer::Type::FixedKla{ kla };
    gi.init_mtr_model(*case_data.simulation, std::move(auto_mtr_type));

    // case_data.simulation->setMtrModel(
    //     Simulation::MassTransfer::Type::FixedKla{kla});

    // case_data.simulation->getter().get_start_time_mut() = start_time;
    case_data.simulation->overwriteStartTime(start_time);
    gi.set_initial_number_particle(np);
    std::cout << "SIMULATION: " << case_data.exec_info.run_id << " LOADED"
              << std::endl;
    return true;
  }

} // namespace SerDe

#endif // USE_CEAREAL
