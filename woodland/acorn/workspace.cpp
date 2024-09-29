#include "woodland/acorn/workspace.hpp"
#include "woodland/acorn/openmp.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cassert>

namespace woodland {
namespace acorn {

int Workspace::get_tid () const {
  const int tid = acorn::in_parallel_region() ? acorn::get_thread_num() : 0;
  assert(size_t(tid) < bufs.size());
  return tid;
}

Workspace::Workspace () { init_threading(); }

void Workspace::init_threading (int nthr) {
  assert(not acorn::in_parallel_region());
  if (nthr < 0) nthr = acorn::get_max_threads();
  if (size_t(nthr) > bufs.size()) bufs.resize(nthr);
}

void* Workspace::set_cap (size_t nbyte) {
  auto& buf = bufs[get_tid()];
  if (buf.size() < nbyte)
    buf.resize(nbyte);
  return buf.data();
}

size_t Workspace::get_cap () const { return bufs[get_tid()].size(); }

void* Workspace::get () { return bufs[get_tid()].data(); }

} // namespace acorn
} // namespace woodland
