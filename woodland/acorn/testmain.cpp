#include "woodland/acorn/unittest.hpp"
#include "woodland/acorn/vv.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/interaction_integrals.hpp"

#include <string>

struct Command {
  enum Enum : int
    { unittest = 0, conv_test_flat_strip, conv_test_circle, conv_test_cylinder,
      conv_test_ellipse, conv_test_ellipse_cylinder, time_calc_sigma_point,
      study_triquad, invalid };
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
  if (e == "conv_test_flat_strip") return conv_test_flat_strip;
  if (e == "conv_test_circle") return conv_test_circle;
  if (e == "conv_test_cylinder") return conv_test_cylinder;
  if (e == "conv_test_ellipse") return conv_test_ellipse;
  if (e == "conv_test_ellipse_cylinder") return conv_test_ellipse_cylinder;
  if (e == "time_calc_sigma_point") return time_calc_sigma_point;
  if (e == "study_triquad") return study_triquad;
  return invalid;
}

std::string Command::convert (const Command::Enum e) {
  switch (e) {
  case unittest: return "unittest";
  case conv_test_flat_strip: return "conv_test_flat_strip";
  case conv_test_circle: return "conv_test_circle";
  case conv_test_cylinder: return "conv_test_cylinder";
  case conv_test_ellipse: return "conv_test_ellipse";
  case conv_test_ellipse_cylinder: return "conv_test_ellipse_cylinder";
  case time_calc_sigma_point: return "time_calc_sigma_point";
  case study_triquad: return "study_triquad";
  case invalid:
  default: return "invalid";
  }
}

bool Command::is_valid (const Enum e) {
  return e >= unittest && e < invalid;
}

int main (int argc, char** argv) {
  using namespace woodland::acorn;
  using namespace vv::convtest;
  Command::Enum command = Command::unittest;
  if (argc > 1)
    command = Command::convert(argv[1]);
  if (command == Command::invalid) {
    printf("%s <command> args...\n", argv[1]);
    return -1;
  }
  printf("#threads %d\n", get_max_threads());
  switch (command) {
  case Command::unittest: return unittest();
  case Command::conv_test_flat_strip: {
    flatstrip::Config c;
    run(c);
  } break;
  case Command::conv_test_circle: {
    gencyl::Config c(gencyl::Config::p_circle);
    run(c);
  } break;
  case Command::conv_test_cylinder: {
    gencyl::Config c(gencyl::Config::p_cylinder);
    run(c);
  } break;
  case Command::conv_test_ellipse: {
    gencyl::Config c(gencyl::Config::p_ellipse);
    run(c);
  } break;
  case Command::conv_test_ellipse_cylinder: {
    gencyl::Config c(gencyl::Config::p_ellipse_cylinder);
    run(c);
  } break;
  case Command::time_calc_sigma_point: {
    fs3d::time_calc_sigma_point(1L << 21);
  } break;
  case Command::study_triquad: {
    fs3d::study_triquad();
  } break;
  case Command::invalid:
  default: {
    printf("Invalid command: %s\n", argv[1]);
    return -1;
  }
  }
  return 0;
}
