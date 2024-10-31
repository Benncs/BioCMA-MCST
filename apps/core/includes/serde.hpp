
#ifdef USE_CEAREAL

#  ifndef __CORE__SERDE_HPP__
#    define __CORE__SERDE_HPP__

#    include <core/case_data.hpp>
#    include <string_view>

namespace SerDe
{
  void save_simulation(const Core::CaseData &case_data);
  bool load_simulation(Core::CaseData &case_data, std::string_view ser_filename);
} // namespace SerDe

#  endif

#endif // USE_CEAREAL