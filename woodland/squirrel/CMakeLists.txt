set(WOODLAND_SQUIRREL_SRCS
  bezier_cubic.cpp
  calc_stress.cpp
  ctzx.cpp
  ctzx_eval.cpp
  ctzx_exact.cpp
  ctzx_surface.cpp
  ctzx_triangulation.cpp
  discretization.cpp
  exact.cpp
  gallery.cpp
  global_z_surface.cpp
  greens_function.cpp
  mesh_based_surface.cpp
  mesh.cpp
  pywrite.cpp
  solver1d.cpp
  testmain.cpp
  unittest.cpp)

set(WOODLAND_SQUIRREL_HEADERS
  bezier_cubic.hpp
  calc_stress.hpp
  ctzx.hpp
  ctzx_surface.hpp
  ctzx_triangulation.hpp
  discretization.hpp
  exact.hpp
  gallery.hpp
  global_z_surface.hpp
  greens_function.hpp
  mesh_based_surface.hpp
  mesh.hpp
  pywrite.hpp
  solver1d.hpp
  squirrel.hpp
  unittest.hpp)

add_library(woodland_squirrel ${WOODLAND_SQUIRREL_SRCS})
if (OpenMP_Fortran_FOUND)
  target_link_libraries(woodland_squirrel OpenMP::OpenMP_Fortran)
endif()
target_include_directories(woodland_squirrel PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>  # for acorn_config.h
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>  # for headers at build time
  $<INSTALL_INTERFACE:include>)                         # for headers from install

add_executable(woodland_squirrel_unittest testmain.cpp)
target_link_libraries(woodland_squirrel_unittest PRIVATE woodland_squirrel woodland_acorn)
add_test(NAME squirrel_unittest COMMAND $<TARGET_FILE:woodland_squirrel_unittest>)

install(TARGETS woodland_squirrel DESTINATION lib EXPORT woodland)
install(TARGETS woodland_squirrel_unittest DESTINATION bin EXPORT woodland)
install(FILES ${WOODLAND_SQUIRREL_HEADERS} DESTINATION include/woodland/squirrel)
