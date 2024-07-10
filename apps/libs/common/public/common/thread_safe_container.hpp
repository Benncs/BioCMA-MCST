#ifndef __COMMON_THREAD_SAFE_CONTAINER_HPP__
#define __COMMON_THREAD_SAFE_CONTAINER_HPP__

#include <cstddef>
#include <vector>

#define USE_OMP_EXECUTOR 1
#define USE_STL_EXEUCTOR 1

#if defined(USE_OMP_EXECUTOR)
#  define __ATOM_INCR__(__val__) _Pragma("omp atomic") __val__++;
#elif defined(USE_STL_EXECUTOR)
// Define behavior for STL executor (e.g., increment without atomic operation)
#  define __ATOM_2INCR__(__val__) __val__++;
#else
#  error "Error: Neither USE_OMP_EXECUTOR nor USE_STL_EXECUTOR is defined"
#endif

#if defined(USE_OMP_EXECUTOR)
#  define __ATOM_DECR__(__val__) _Pragma("omp atomic") __val__--;
#elif defined(USE_STL_EXECUTOR)
// Define behavior for STL executor (e.g., increment without atomic operation)
#  define __ATOM_2INCR__(__val__) __val__++;
#else
#  error "Error: Neither USE_OMP_EXECUTOR nor USE_STL_EXECUTOR is defined"
#endif



// template<typename T,typename F_merge>
// class ThreadSafeDataContainer
// {
//     public:

//     T& operator[](size_t );
//     T operator[](size_t )const;
//     ~ThreadSafeDataContainer();
//     ThreadSafeDataContainer(const ThreadSafeDataContainer&)=delete;
//     ThreadSafeDataContainer(ThreadSafeDataContainer&&)=delete;
//     ThreadSafeDataContainer& operator=(ThreadSafeDataContainer&&)=delete;
//     ThreadSafeDataContainer& operator=(const
//     ThreadSafeDataContainer&)=delete;

//     private:
//     // T _data;
//     std::vector<T> _data;
//     F_merge f_ptr;
// };

#endif //__COMMON_THREAD_SAFE_CONTAINER_HPP__