cmake_minimum_required(VERSION 3.21)

project(
  woodland
  VERSION 0.0
  LANGUAGES CXX Fortran)

include(CTest)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
  link_libraries(OpenMP::OpenMP_CXX)
endif()
if (OpenMP_Fortran_FOUND)
  link_libraries(OpenMP::OpenMP_Fortran)
endif()

add_subdirectory(woodland/acorn)

install(EXPORT woodland DESTINATION cmake)
install(FILES cmake/woodland-config.cmake DESTINATION cmake)
