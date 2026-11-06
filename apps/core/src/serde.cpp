

#include "mixture/species_descriptor.hpp"
#include "simulation/descriptors/dimensions.hpp"
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

namespace
{
  std::string
  format_serde_file_path(std::string_view path, std::size_t rank)
  {

    return IO::format(path, std::to_string(rank), ".raw");
  }

  void
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

  void
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

  struct SerDePayload
  {

    std::string version;
    ExecInfo serde_exec{};
    Mixture::SpecieTable t;
    uint64_t np = 0;
    Simulation::Dimensions dims;
    std::vector<double> read_c_liq;
    std::optional<std::vector<double>> read_c_gas;
    double start_time{};
    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    SerDePayload(ExecInfo exec,
                 Mixture::SpecieTable table,
                 uint64_t n_p,
                 Simulation::Dimensions dimension,
                 std::vector<double> c_liq,
                 std::optional<std::vector<double>> c_gas,
                 double stime,
                 std::unique_ptr<MC::MonteCarloUnit> mcu)
        : serde_exec(exec), t(std::move(table)), np(n_p), dims(dimension),
          read_c_liq(std::move(c_liq)), read_c_gas(std::move(c_gas)),
          start_time(stime), mc_unit(std::move(mcu))
    {
      this->version = ExecInfo::get_version();
    }

    SerDePayload() = default;

    template <typename Archive>
    void
    load(Archive& ar)
    {
      ar(version,
         serde_exec,
         t,
         np,
         dims,
         read_c_liq,
         read_c_gas,
         start_time,
         mc_unit);
    }

    template <typename Archive>
    void
    save(Archive& ar) const
    {

      ar(cereal::make_nvp("version", version));
      ar(serde_exec, t, np, dims, read_c_liq, read_c_gas, start_time, mc_unit);
    }

    template <typename Archive>
    static SerDePayload
    get(Archive& ar)
    {
      SerDePayload payload;
      ar(payload);
      return payload;
    }
  };

} // namespace

namespace SerDe
{

  // using Archive_t = cereal::XMLOutputArchive;
  //  using iArchive_t = cereal::XMLInputArchive;
  using Archive_t = cereal::BinaryOutputArchive;
  using iArchive_t = cereal::BinaryInputArchive;
  void
  save_simulation(const Core::CaseData& case_data)
  {

    const auto serde_name = IO::format(case_data.params.results_file_name,
                                       "_serde_",
                                       case_data.exec_info.current_rank,
                                       ".raw");

    auto accessor = case_data.simulation->getter();
    auto dim = accessor.getDimensions();
    auto cliq = accessor.getCliqData();
    auto t = accessor.table_ownable();
    auto cgas = accessor.getCgasData();
    auto mc_unit = case_data.simulation->move_mcu();

    auto cgas_a = cgas.has_value() ? std::make_optional(std::vector<double>(
                                         cgas->begin(), cgas->end()))
                                   : std::nullopt;

    std::vector<double> c_liq_save(cliq.begin(), cliq.end());

    SerDePayload payload(case_data.exec_info,
                         std::move(*t),
                         case_data.params.number_particle,
                         dim,
                         std::move(c_liq_save),
                         std::move(cgas_a),
                         accessor.absolute_time(),
                         std::move(mc_unit));

    std::ostringstream buf(std::ios::binary);
    Archive_t ar(buf);
    ar(payload);
    write_to_file(buf, serde_name);
  }

  namespace
  {

    Simulation::ScalarInitializer
    build_scalar_init(Core::GlobalInitialiser& gi,
                      const Simulation::Dimensions& dims,
                      std::vector<double>&& c_liq_buffer,
                      std::optional<std::vector<double>>&& c_gas_buffer)
    {

      // TODO: Scalar init is not optional
      auto sc = gi.init_scalar();

      /*
       *
       * At this stage, init scalar should not be optional
       */

      if (!sc.has_value())
      {

        throw std::invalid_argument("Scalar init is not optional");
      }

      if (c_liq_buffer.size() % sc->volumesliq.size() != 0)
      {
        throw std::invalid_argument("Liquid Concentration buffer and CM volume "
                                    "have incompatible sizes ");
      }

      if (c_gas_buffer)
      {
        if (c_gas_buffer->size() % sc->volumesgas.size() != 0)
        {
          throw std::invalid_argument("Gas Concentration buffer and CM volume "
                                      "have incompatible sizes ");
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

      return *sc;
    }

    bool
    load_simulation(Core::GlobalInitialiser& gi,
                    Core::CaseData& case_data,
                    std::string_view ser_filename)
    {

      std::stringstream buffer;
      read_file(buffer, ser_filename);

      iArchive_t ar(buffer);

      auto [version,
            serde_exec,
            t,
            np,
            dims,
            read_c_liq,
            read_c_gas,
            start_time,
            mc_unit]
          = SerDePayload::get(ar);

      assert(mc_unit != nullptr);

      if (mc_unit == nullptr)
      {

        std::cout << "MCU load failed" << std::endl;
        return false;
      }

      case_data.exec_info.run_id = serde_exec.run_id;

      auto sc = build_scalar_init(
          gi, dims, std::move(read_c_liq), std::move(read_c_gas));

      gi.set_initial_number_particle(np);
      gi.set_table(std::make_shared<Mixture::SpecieTable>(t));
      auto simulation = gi.init_simulation(std::move(mc_unit), std::move(sc));

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
      case_data.simulation->overwriteStartTime(start_time);

      std::cout << "SIMULATION: " << case_data.exec_info.run_id << " LOADED"
                << std::endl;
      return true;
    }

  } // namespace

  std::optional<Core::CaseData>
  impl_load(
      std::shared_ptr<IO::Logger>& logger,
      const ExecInfo& exec,
      const Core::UserControlParameters&& params,
      std::optional<Simulation::Feed::SimulationFeed> feed,
      std::optional<Simulation::MassTransfer::Type::MtrTypeVariant> mtr_type)
  {
    Core::CaseData case_data;
    case_data.exec_info = exec;

    Core::GlobalInitialiser gi(exec, params, nullptr, logger);
    auto transition = gi.init_transitionner();

    if (!transition.has_value())
    {
      return std::nullopt;
    }

    case_data.transitioner = std::move(*transition);
    gi.init_feed(std::move(feed));
    gi.init_mtr_model(std::move(mtr_type));

    if (!params.serde_file.has_value())
    {
      return std::nullopt;
    }

    const auto serde_filename
        = format_serde_file_path(*params.serde_file, exec.current_rank);

    if (logger)
    {
      logger->print("Serde", serde_filename);
    }

    try
    {
      const bool ok_init
          = SerDe::load_simulation(gi, case_data, serde_filename);

      if (!gi.check_init_terminate() || !ok_init)
      {
        return std::nullopt;
      }
    }
    catch (std::exception& e)
    {
      auto err = "CORE::load::load_simulation:" + std::string(e.what());
      throw std::runtime_error(err);
    }

    case_data.params = gi.get_parameters();
    // TODO: Check integreity between transitionner and read value
    // (transitionner and simulation tpf
    // + n_species)

    return case_data;
  }

  std::unique_ptr<MC::MonteCarloUnit>
  load_mc(std::string_view path)
  {
    std::stringstream buffer;
    read_file(buffer, path);
    iArchive_t ar(buffer);

    auto payload = SerDePayload::get(ar);
    return std::move(payload.mc_unit);
  }

  void
  swap_mc(std::unique_ptr<MC::MonteCarloUnit>&& mcu,
          std::string_view ref_path,
          std::string_view new_path)
  {
    std::stringstream buffer;
    read_file(buffer, ref_path);
    iArchive_t ar(buffer);

    auto payload = SerDePayload::get(ar);
    payload.mc_unit = std::move(mcu);

    std::ostringstream buf(std::ios::binary);
    Archive_t ar_w(buf);
    ar_w(payload);
    write_to_file(buf, new_path);
  }

} // namespace SerDe

#endif // USE_CEAREAL
