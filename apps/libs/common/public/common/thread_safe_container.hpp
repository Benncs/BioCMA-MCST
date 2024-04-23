#ifndef __COMMON_THREAD_SAFE_CONTAINER_HPP__
#define __COMMON_THREAD_SAFE_CONTAINER_HPP__

#include <vector>
#include <cstddef>


#define __ATOM_INCR__(__val__) {\
_Pragma("omp atomic") \
__val__++;\
}

#define __ATOM_DECR__(__val__) {\
_Pragma("omp atomic") \
__val__--;\
}




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
//     ThreadSafeDataContainer& operator=(const ThreadSafeDataContainer&)=delete;

//     private:
//     // T _data;
//     std::vector<T> _data;
//     F_merge f_ptr;
// };



#endif //__COMMON_THREAD_SAFE_CONTAINER_HPP__