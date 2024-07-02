#include <model_list.hpp>
#include <string>
#include <unordered_map>
#include <stdexcept>

#ifdef USE_PYTHON_MODULE
#  include <pymodule/import_py.hpp>
#endif

@INCLUDES@


KModel load_model_(const std::string& name, bool use_python_module)
{
   #ifdef USE_PYTHON_MODULE
    if (use_python_module)
    {
        std::string module_name = std::string("modules.") + name;
        return get_python_module(module_name);
     
    }
    #endif
    @BODY@
}

std::vector<std::string> get_available_models(bool use_python_module)
{
    std::vector<std::string> list;
    @AM_BODY@

    return list;
}