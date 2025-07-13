#ifndef INCLUDE_WOODLAND_ACORN_MACROS
#define INCLUDE_WOODLAND_ACORN_MACROS

#include <sstream>
#include <stdexcept>

#include "woodland/acorn/dbg.hpp"

// The following macros print to stderr as well as throw std::logic_error. To
// prevent the print to stderr, set acorn::print_throw_if_msg_to_stderr = false.

namespace woodland {
namespace acorn {
extern bool print_throw_if_msg_to_stderr;
}
}

#define throw_if_nomsg(condition) do {                                  \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\n is an error";                              \
      if (woodland::acorn::print_throw_if_msg_to_stderr)                \
        fprintf(stderr, "%s\n", _ss_.str().c_str());                    \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

#define throw_if(condition, message) do {                               \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
      if (woodland::acorn::print_throw_if_msg_to_stderr)                \
        fprintf(stderr, "%s\n", _ss_.str().c_str());                    \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

// Wrapper to assert for routines that return a boolean, with true being OK and
// false an unexpected result. Using this, we don't need to put NDEBUG code at
// the call site to prevent unused-variable warnings when NDEBUG is defined.
#ifdef NDEBUG
# define assert_ok(ok) ok
#else
# define assert_ok(ok) assert(ok)
#endif

#define rununittest(f) do {                     \
    const auto t0 = acorn::dbg::gettime();      \
    ne = f();                                   \
    const auto t1 = acorn::dbg::gettime();      \
    if (ne) printf(#f " FAILED: ne %d\n", ne);  \
    printf("%60s: %6.3f\n", #f, t1 - t0);       \
    fflush(stdout);                             \
    nerr += ne;                                 \
  } while (0)

#endif
