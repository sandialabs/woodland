set(WOODLAND_ACORN_SRCS
  compose_quadrature.cpp
  compose_triquad.cpp
  elastostatics.cpp
  elastostatics_integrals.cpp
  fs3d.cpp
  fs3d_sigma_znml.cpp
  hfp.cpp
  hs3d.cpp
  hs3d_sigma_A.cpp
  hs3d_sigma_B.cpp
  hs3d_sigma_C.cpp
  interaction_integrals.cpp
  linalg.cpp
  matvec.cpp
  plane_geometry.cpp
  triangle.cpp
  unittest.cpp
  util.cpp
  vv.cpp
  vv_flat_strip.cpp
  workspace.cpp)

set(WOODLAND_ACORN_HEADERS
  acorn.hpp
  caller_integrand.hpp
  compose_lagrange_polynomial.hpp
  compose_quadrature.hpp
  compose_triquad.hpp
  dbg.hpp
  elastostatics.hpp
  elastostatics_integrals.hpp
  fs3d.hpp
  hfp.hpp
  hs3d.hpp
  interaction_integrals.hpp
  linalg.hpp
  macros.hpp
  matvec.hpp
  openmp.hpp
  plane_geometry.hpp
  triangle.hpp
  types.hpp
  unittest.hpp
  util.hpp
  vectorization.hpp
  vv.hpp
  workspace.hpp)

# Not currently supported:
#   option(WOODLAND_ACORN_VECTORIZE "Enable vectorization." OFF)
set(WOODLAND_ACORN_PACK_N 4 CACHE STRING "Pack size for vectorization.")

set(WOODLAND_ACORN_DC3D_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../extern/dc3d.f")
if (EXISTS ${WOODLAND_ACORN_DC3D_SRC})
  file(READ ${WOODLAND_ACORN_DC3D_SRC} WOODLAND_ACORN_DC3D_TEXT)
  string(FIND "${WOODLAND_ACORN_DC3D_TEXT}" "CHARACTER SPACE" WOODLAND_ACORN_OUTPUT)
  if (WOODLAND_ACORN_OUTPUT EQUAL -1)
    message("acorn> extern/dc3d.f was found, but the patch has not been applied.\n \
      Several unit tests are not available.\n \
      See extern/README.md for instructions.")
  else()
    list(APPEND WOODLAND_ACORN_SRCS ${WOODLAND_ACORN_DC3D_SRC})
    set(WOODLAND_ACORN_HAVE_DC3D TRUE)
    message("acorn> extern/dc3d.f was found, and the patch was applied successfully.\n \
      All unit tests are available.")
  endif()
else()
  message("acorn> extern/dc3d.f was not found; building without it.\n \
      See extern/README.md for instructions to obtain it.\n \
      Several unit tests are not available.")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acorn_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/acorn_config.h)

add_library(woodland_acorn ${WOODLAND_ACORN_SRCS})
if (OpenMP_Fortran_FOUND)
  target_link_libraries(woodland_acorn OpenMP::OpenMP_Fortran)
endif()
target_include_directories(woodland_acorn PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>  # for acorn_config.h
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>  # for headers at build time
  $<INSTALL_INTERFACE:include>)                         # for headers from install

add_executable(woodland_acorn_unittest testmain.cpp)
target_link_libraries(woodland_acorn_unittest PRIVATE woodland_acorn)
add_test(NAME acorn_unittest COMMAND $<TARGET_FILE:woodland_acorn_unittest>)

install(TARGETS woodland_acorn DESTINATION lib EXPORT woodland)
install(TARGETS woodland_acorn_unittest DESTINATION bin EXPORT woodland)
install(FILES ${WOODLAND_ACORN_HEADERS} DESTINATION include/woodland/acorn)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/acorn_config.h DESTINATION include/woodland/acorn)
