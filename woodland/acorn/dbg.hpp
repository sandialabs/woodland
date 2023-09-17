#ifndef INCLUDE_WOODLAND_ACORN_DBG
#define INCLUDE_WOODLAND_ACORN_DBG

#include <sys/time.h>

#include <cstdlib>
#include <iostream>
#include <sstream>

namespace woodland {
namespace acorn {

#define forl(k,klim) for (int k = 0; k < klim; ++k)

namespace dbg {
static inline int getrank () {
  int pid = 0;
#if ! defined WOODLAND_ACORN_NO_MPI && defined _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
  return pid;
}
} // namespace dbg

#ifndef pr
#define pr(m) do {                                    \
    int _pid_ = acorn::dbg::getrank();                \
    std::stringstream _ss_;                           \
    _ss_ << "pid " << _pid_ << " " << m << std::endl; \
    std::cerr << _ss_.str();                          \
  } while (0)
#define pr0(m) do {                                     \
    int _pid_ = acorn::dbg::getrank();                  \
    MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);              \
    if (_pid_ != 0) break;                              \
    std::stringstream _ss_;                             \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define pr0c(m) pr0(#m << " | " << (m))
#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
#endif

template <typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::stringstream ss;
  ss << name << " = [";
  for (size_t i = 0; i < n; ++i) ss << " " << v[i];
  ss << "]";
  pr(ss.str());
}
template <typename Array>
static void prarr (const std::string& name, const Array& a) {
  prarr(name, a.data(), a.size());
}

#ifndef mprarr
# define mprarr(a) prarr(#a, a)
#endif

namespace dbg {
template <typename T> T strto(const char* s);
template <> inline int strto (const char* s) { return std::atoi(s); }
template <> inline bool strto (const char* s) { return std::atoi(s); }
template <> inline double strto (const char* s) { return std::atof(s); }
template <> inline std::string strto (const char* s) { return std::string(s); }

template <typename T>
bool getenv (const std::string& varname, T& var) {
  const char* var_s = std::getenv(varname.c_str());
  if ( ! var_s) return false;
  var = strto<T>(var_s);
  return true;
}
} // namespace dbg

namespace dbg {
static inline double gettime () {
#if 0
  return 0;
#else
  timeval t;
  gettimeofday(&t, 0);
  static const double us = 1e6;
  return (t.tv_sec*us + t.tv_usec)/us;
#endif
}
} // namespace dbg

#define sumtime(name, t0, t1) (name) += (t1) - (t0)
//#define sumtime(name, t0, t1)

} // namespace acorn
} // namespace woodland

#endif
