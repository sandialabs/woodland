#ifndef INCLUDE_WOODLAND_ACORN_MACROS
#define INCLUDE_WOODLAND_ACORN_MACROS

#include <sstream>
#include <stdexcept>

#define throw_if_nomsg(condition) do {                                  \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\n is an error";                              \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

#define throw_if(condition, message) do {                               \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
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

#endif
