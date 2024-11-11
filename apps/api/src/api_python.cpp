#include <pybind11/pybind11.h>
#include <api/api_raw.h>  // Include the header with your C functions.
#include <api/api.hpp>
#include <stdexcept>
namespace py = pybind11;

PYBIND11_MODULE(handle_module, m) {
    // Wrapping the Handle structure
    py::class_<Handle>(m, "Handle");
    m.def("init_handle_raw", [](uint32_t n_rank,uint32_t i_rank){
        auto opt = Handle::init(n_rank, i_rank);
        if(opt.has_value())
        {
            return std::move(opt.value());
        }
            throw std::runtime_error("Simulation handle initialisation failed");
    });

    m.def("load", [](uint32_t n_rank,uint32_t i_rank){
        auto opt = Handle::load(n_rank, i_rank);
        if(opt.has_value())
        {
            return std::move(opt.value());
        }
            throw std::runtime_error("Simulation handle initialisation failed");
    });

    m.def("exec", &exec);
    m.def("register_parameters", &register_parameters);
    m.def("apply", &apply);

}
