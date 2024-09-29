#ifndef INCLUDE_WOODLAND_ACORN_WORKSPACE
#define INCLUDE_WOODLAND_ACORN_WORKSPACE

#include <memory>
#include <vector>

namespace woodland {
namespace acorn {

struct Workspace {
  typedef std::shared_ptr<Workspace> Ptr;
  
  Workspace();

  // Initialize outside of a threaded region. If nthr is not provide,
  // get_max_threads() is used.
  void init_threading(int nthr = -1);

  // These work inside or outside a parallel region.
  void* set_cap(size_t nbyte);
  size_t get_cap() const;
  void* get();

private:
  std::vector<std::vector<char>> bufs;

  inline int get_tid() const;
};

} // namespace acorn
} // namespace woodland

#endif
