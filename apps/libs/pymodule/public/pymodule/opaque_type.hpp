#ifndef __PY_OPAQUE_HPP__
#define __PY_OPAQUE_HPP__

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
struct OpaquePointer
{
  void *ptr{};
  OpaquePointer() = default;
  OpaquePointer(const OpaquePointer &) = delete;
  OpaquePointer(OpaquePointer &&rhs) noexcept : ptr(rhs.ptr)
  {

    rhs.ptr = nullptr;
  };

  OpaquePointer &operator=(OpaquePointer &&rhs)
 noexcept   {
    if (&rhs != this)
    {
      ptr = rhs.ptr;
      rhs.ptr = nullptr;
    }
    return *this;
  };
  OpaquePointer &operator=(const OpaquePointer &&) = delete;

  ~OpaquePointer()
  {
    if (ptr != nullptr)
    {
      auto *obj = static_cast<py::object *>(ptr);
      delete obj;
    }
  }
};

void declare_opaque(py::module &m);



#endif //__PY_OPAQUE_HPP__
