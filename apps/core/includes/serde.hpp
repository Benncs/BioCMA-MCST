#ifdef USE_CEAREAL

#ifndef __CORE__SERDE_HPP__
#define __CORE__SERDE_HPP__


#include <core/case_data.hpp>
namespace SerDe
{
    void save_simulation(const Core::CaseData& case_data);
} //namespace SerDe

#endif 

#endif // USE_CEAREAL