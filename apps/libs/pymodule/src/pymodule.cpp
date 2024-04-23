#include <any>
#include <iostream>
#include <mc/particles/mcparticles.hpp>
#include <memory>
#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pymodule/opaque_type.hpp>
namespace py = pybind11;

PYBIND11_MODULE(pyBioCMAMCST, m)
{

  py::enum_<MC::CellStatus>(m, "CellStatus")
      .value("IDLE", MC::CellStatus::IDLE)
      .value("DEAD", MC::CellStatus::DEAD)
      .value("CYTOKINESIS", MC::CellStatus::CYTOKINESIS);

  py::class_<MC::Particles>(m, "Particles")
      .def_readwrite("current_container",
                     &MC::Particles::current_container,
                     py::return_value_policy::reference)
      .def_readwrite("current_domain", &MC::Particles::current_domain)
      .def_readwrite("random_seed", &MC::Particles::random_seed)
      .def_readwrite("id", &MC::Particles::id)
      .def_readwrite("weight", &MC::Particles::weight)
      .def_readwrite("status", &MC::Particles::status)
      // .def("child",
      //      [](MC::Particles &p)
      //      {
      //        MC::Particles child(p);
      //        return child;
      //      })

      .def("getOpaque",
           [](MC::Particles &p)
           { return std::any_cast<std::shared_ptr<OpaquePointer>>(p.data); });

  declare_opaque(m);
}