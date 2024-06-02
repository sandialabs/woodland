set(WOODLAND_CONVZX_SRCS
  calc_stress.cpp
  convtest_zx.cpp
  convtest_zx_eval.cpp
  convtest_zx_exact.cpp
  convtest_zx_surface.cpp
  convtest_zx_hmmvp.cpp
  discretization.cpp
  exact.cpp
  gallery.cpp
  pywrite.cpp
  testmain.cpp
  triangulation.cpp
  unittest.cpp)

add_executable(woodland_examples_convzx ${WOODLAND_CONVZX_SRCS})
target_link_libraries(woodland_examples_convzx PRIVATE woodland_acorn)
if (OpenMP_Fortran_FOUND)
  target_link_libraries(woodland_examples_convzx PRIVATE OpenMP::OpenMP_Fortran)
endif()
add_test(NAME examples_convzx_unittest COMMAND $<TARGET_FILE:woodland_examples_convzx>)
install(TARGETS woodland_examples_convzx DESTINATION bin EXPORT woodland)
