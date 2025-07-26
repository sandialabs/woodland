#include "woodland/squirrel/greens_function.hpp"

namespace woodland {
namespace squirrel {

bool operator== (const GreensFnParams& a, const GreensFnParams& b) {
  return (a.lam == b.lam and a.mu == b.mu and a.halfspace == b.halfspace);
}

bool operator!= (const GreensFnParams& a, const GreensFnParams& b) {
  return not (a == b);
}

} // namespace squirrel
} // namespace woodland
