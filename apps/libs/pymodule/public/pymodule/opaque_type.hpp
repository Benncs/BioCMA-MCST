#ifndef __PY_OPAQUE_HPP__
#define __PY_OPAQUE_HPP__

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
struct OpaquePointer
{
  void* ptr{nullptr}; // Ensure ptr is initialized to nullptr by default

  OpaquePointer() = default; // Default constructor

  OpaquePointer(const OpaquePointer& rhs) noexcept =
      default; // Default copy constructor

  OpaquePointer& operator=(const OpaquePointer& rhs) = default;

  OpaquePointer(OpaquePointer&& rhs) noexcept : ptr(rhs.ptr)
  { // Move constructor
    rhs.ptr = nullptr;
  }

  OpaquePointer& operator=(OpaquePointer&& rhs) noexcept
  { // Move assignment operator
    if (this != &rhs)
    { // Prevent self-assignment
      ptr = rhs.ptr;
      rhs.ptr = nullptr;
    }
    return *this;
  }

  // Add a destructor to handle cleanup if necessary
  ~OpaquePointer() = default;
};

void declare_opaque(py::module& m);

#endif //__PY_OPAQUE_HPP__
