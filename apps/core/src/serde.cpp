#ifdef USE_CEAREAL
#include "common/execinfo.hpp"
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/array.hpp> //MC::events use vector internally
#include <cereal/types/optional.hpp>
#include <cereal/types/string.hpp> //MC::List use vector internally
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp> //MC::List use vector internally
#include <fstream>
#include <optional>
#include <serde.hpp>
#include <string_view>

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

#define PLAINTEXT_ARCHIVE cereal::XMLOutputArchive
#define BINARY_ARCHIVE cereal::BinaryArchive

using Archive_t = cereal::XMLOutputArchive;

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
      auto dim = case_data.simulation->getDim();
      auto cliq = case_data.simulation->getCliqData();
      auto cgas = case_data.simulation->getCgasData();

      std::optional<std::vector<double>> cgas_a =
          cgas.has_value() ? std::make_optional(std::vector<double>(cgas->begin(), cgas->end())) : std::nullopt;

      ar(dim, std::vector<double>(cliq.begin(), cliq.end()), cgas_a);

      ar(*case_data.simulation->mc_unit);
    }

    write_to_file(buf, serde_name.str());
  }
} // namespace SerDe


#endif // USE_CEAREAL