#ifndef __EIGEN_DIAG_HPP__
#define __EIGEN_DIAG_HPP__

#define DO_PRAGMA(x) _Pragma(#x)

#ifndef IGNORE_EIGEN_DIAG

#  define EIGEN_DIAG_PUSH                                                      \
    DO_PRAGMA(GCC diagnostic push)                                             \
    DO_PRAGMA(GCC diagnostic ignored "-Wunused-but-set-variable")              \
    DO_PRAGMA(GCC diagnostic ignored "-Wnan-infinity-disabled")

#  define EIGEN_DIAG_POP DO_PRAGMA(GCC diagnostic pop)

#else

/* Diagnostics disabled → macros expand to nothing */
#  define EIGEN_DIAG_PUSH
#  define EIGEN_DIAG_POP

#endif

#endif