#ifndef INCLUDE_WOODLAND_ACORN_OPENMP
#define INCLUDE_WOODLAND_ACORN_OPENMP

#ifdef _OPENMP
# include <omp.h>
# define ompparfor  _Pragma("omp parallel for")
# define omppar     _Pragma("omp parallel")
# define ompfor     _Pragma("omp for")
# define ompbarrier _Pragma("omp barrier")
# define ompsingle  _Pragma("omp single")
# define ompmaster  _Pragma("omp master")
# define ompatomic  _Pragma("omp atomic")
# define ompsimd    _Pragma("omp simd")
#else
# define ompparfor
# define omppar
# define ompfor
# define ompbarrier
# define ompsingle
# define ompmaster
# define ompatomic
# define ompsimd
#endif

namespace woodland {
namespace acorn {

inline int get_max_threads () {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline int get_thread_num () {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif  
}

inline bool in_parallel_region () {
#ifdef _OPENMP
  return omp_in_parallel();
#else
  return false;
#endif
}

} // namespace acorn
} // namespace woodland

#endif
