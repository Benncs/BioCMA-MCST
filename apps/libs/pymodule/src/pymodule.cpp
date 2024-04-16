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

  py::class_<MC::Particles>(m, "Particles")
      .def_readwrite("current_container",
                     &MC::Particles::current_container,
                     py::return_value_policy::reference)
      .def_readwrite("current_domain", &MC::Particles::current_domain)
      .def_readwrite("random_seed", &MC::Particles::random_seed)
      .def_readwrite("id", &MC::Particles::id)
      .def_readwrite("weight", &MC::Particles::weight)

      .def("getOpaque",
           [](MC::Particles &p)
           { return std::any_cast<std::shared_ptr<OpaquePointer>>(p.data); });
  // .def("setmodel",
  //      [](MC::Particles &p, const PythonCustomModel &new_model)
  //      { p.data = new_model; })
  // .def("getModel",
  //      [](MC::Particles &p)
  //      { return std::any_cast<PythonCustomModel>(p.data); });

  py::class_<OpaquePointer, std::shared_ptr<OpaquePointer>>(
      m, "OpaquePointer", py::call_guard<py::gil_scoped_release>())
      .def(py::init<>())

      .def_readwrite("ptr", &OpaquePointer::ptr)
      .def(
          "init",
          [](OpaquePointer &opaque_ptr, py::object obj)
          { opaque_ptr.ptr = new py::object(obj); },
          py::return_value_policy::reference_internal)

      .def(
          "cast",
          [](OpaquePointer &opaque_ptr)
          {
            if (opaque_ptr.ptr)
            {
              py::object *obj_ptr = static_cast<py::object *>(opaque_ptr.ptr);
              return py::reinterpret_borrow<py::object>(*obj_ptr);
            }
            else
            {
              return py::reinterpret_borrow<py::object>(nullptr);
            }
          },
          pybind11::return_value_policy::reference);
}