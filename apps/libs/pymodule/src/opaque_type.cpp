#include <pymodule/opaque_type.hpp>

void declare_opaque(py::module& m)
{
  py::class_<OpaquePointer, std::shared_ptr<OpaquePointer>>(
      m, "OpaquePointer", py::call_guard<py::gil_scoped_release>())
      .def(py::init<>())

      // .def_readwrite("ptr", &OpaquePointer::ptr)
      // .def(
      //     "init",
      //     [](OpaquePointer &opaque_ptr, py::object obj)
      //     { opaque_ptr.ptr = new py::object(obj); },
      //     py::return_value_policy::reference_internal)
      .def("__repr__",
           [](OpaquePointer& opaque_ptr) -> std::string
           {
             if (opaque_ptr.ptr != nullptr)
             {
               return "<OpaquePointer holding Python object>";
             }
             return "<OpaquePointer empty>";
           })

      .def(
          "get",
          [](OpaquePointer& opaque_ptr)
          {
            if (opaque_ptr.ptr != nullptr)
            {
              auto* obj_ptr = static_cast<py::object*>(opaque_ptr.ptr);
              return py::reinterpret_borrow<py::object>(*obj_ptr);
            }

             return py::reinterpret_borrow<py::object>(py::none());
          },
          pybind11::return_value_policy::reference);
}