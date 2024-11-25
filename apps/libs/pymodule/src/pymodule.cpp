#include "mc/particles/data_holder.hpp"
#include "pymodule/opaque_type.hpp"
#include <cstdlib>
#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/gil.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pymodule/import_py.hpp>
#include <pymodule/opaque_type.hpp>

PYBIND11_EMBEDDED_MODULE(cpp_module, m) {
    m.attr("a") = 1;
    declare_opaque(m);
    pybind11::class_<MC::ParticleDataHolder>(m,"ParticleDataHolder")
    .def_readwrite("id", &MC::ParticleDataHolder::id);

}


// PYBIND11_MODULE(pyBioCMAMCST, m)
// {

//   py::enum_<MC::CellStatus>(m, "CellStatus")
//       .value("IDLE", MC::CellStatus::IDLE)
//       .value("DEAD", MC::CellStatus::DEAD)
//       .value("CYTOKINESIS", MC::CellStatus::CYTOKINESIS);

//   py::class_<MC::Particles>(m, "Particles")
//       .def_readwrite("current_container",
//                      &MC::Particles::current_container,
//                      py::return_value_policy::reference)
//       .def_readwrite("current_domain", &MC::Particles::current_domain)
//       .def_readwrite("random_seed", &MC::Particles::random_seed)
//       .def_readwrite("id", &MC::Particles::id)
//       .def_readwrite("weight", &MC::Particles::weight)
//       .def_readwrite("status", &MC::Particles::status)
//       // .def("child",
//       //      [](MC::Particles &p)
//       //      {
//       //        MC::Particles child(p);
//       //        return child;
//       //      })

//       .def("getOpaque",
//            [](MC::Particles &p)
//            { return std::any_cast<std::shared_ptr<OpaquePointer>>(p.data); });

//   declare_opaque(m);
// }