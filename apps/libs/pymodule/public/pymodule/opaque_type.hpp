#ifndef __PY_OPAQUE_HPP__
#define __PY_OPAQUE_HPP__

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
struct OpaquePointer
{
  void *ptr;

  ~OpaquePointer()
  {
    if(ptr)
    {
      py::object* obj = static_cast<py::object*>(ptr);
      delete obj;
    }
  }
};


#define GIL_RELEASE pybind11::gil_scoped_release release;
#define GIL_ACQUIRE pybind11::gil_scoped_acquire acquire;

#define PY_SAFE_CONTEXT(expr)                                                  \
  {                                                                            \
    Py_BEGIN_ALLOW_THREADS ;gstate = PyGILState_Ensure();                       \
    expr PyGILState_Release(gstate);                                           \
    Py_END_ALLOW_THREADS                                                       \
  }

#endif //__PY_OPAQUE_HPP__
