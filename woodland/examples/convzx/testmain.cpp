#include "woodland/acorn/openmp.hpp"
#include "woodland/examples/convzx/unittest.hpp"
#include "woodland/examples/convzx/convtest_zx.hpp"

struct Command {
  // ct_ab_zx: convergence test of a w.r.t. reference b, for extruded z(x)
  // surface. w is Woodland, o is Okada.
  enum Enum : int
    { unittest = 0, runcase,
      ct_we_zx, ct_oe_zx,
      invalid };
  static Enum convert(const int i);
  static Enum convert(const std::string& e);
  static std::string convert(const Enum e);
  static bool is_valid(const Enum e);
};

Command::Enum Command::convert (const int i) {
  auto s = static_cast<Command::Enum>(i);
  if ( ! is_valid(s)) s = invalid;
  return s;
}

Command::Enum Command::convert (const std::string& e) {
  if (e == "unittest") return unittest;
  if (e == "runcase") return runcase;
  if (e == "ct_we_zx") return ct_we_zx;
  if (e == "ct_oe_zx") return ct_oe_zx;
  return invalid;
}

std::string Command::convert (const Command::Enum e) {
  switch (e) {
  case unittest: return "unittest";
  case runcase: return "runcase";
  case ct_we_zx: return "ct_we_zx";
  case ct_oe_zx: return "ct_oe_zx";
  case invalid:
  default: return "invalid";
  }
}

bool Command::is_valid (const Enum e) {
  return e >= unittest && e < invalid;
}

int main (int argc, char** argv) {
  using namespace woodland::examples::convzx;
  Command::Enum command = Command::unittest;
  if (argc > 1)
    command = Command::convert(argv[1]);
  if (command == Command::invalid) {
    printf("Invalid command: %s\n", argv[1]);
    printf("%s <command> args...\n", argv[0]);
    return -1;
  }
  printf("#threads %d\n", woodland::acorn::get_max_threads());
  switch (command) {
  case Command::unittest: return unittest();
  case Command::runcase:
  case Command::ct_we_zx:
  case Command::ct_oe_zx: {
    std::string params = "testcase=0,element=spline-nml-4,disloc=2";
    if (argc > 2) params = argv[2];
    switch (command) {
    case Command::runcase: run_case(params); break;
    case Command::ct_we_zx: convtest_w_vs_e(params); break;
    case Command::ct_oe_zx: convtest_o_vs_e(params); break;
    default: assert(0);
    }
  } break;
  case Command::invalid:
  default: {
    printf("Invalid command: %s\n", argv[1]);
    return -1;
  }
  }
  return 0;
}
