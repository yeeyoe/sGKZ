cmake_minimum_required(VERSION 3.20)
project(shortest_gkz LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(gkz_core
  src/gkz.cpp
)
target_include_directories(gkz_core PUBLIC include)
target_link_libraries(gkz_core PUBLIC CGAL::CGAL Eigen3::Eigen)
target_compile_options(gkz_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(shortest_gkz src/main.cpp)
target_link_libraries(shortest_gkz PRIVATE gkz_core)

include(CTest)
if(BUILD_TESTING)
  add_executable(gkz_tests tests/gkz_tests.cpp)
  target_link_libraries(gkz_tests PRIVATE gkz_core)
  add_test(NAME gkz_tests COMMAND gkz_tests)
  add_executable(oracle_kernel_consistency_test
    tests/oracle_kernel_consistency_test.cpp)
  target_link_libraries(oracle_kernel_consistency_test PRIVATE gkz_core)
  add_test(NAME oracle_kernel_consistency
    COMMAND oracle_kernel_consistency_test)
  add_test(
    NAME cli_square
    COMMAND shortest_gkz
      --points ${CMAKE_CURRENT_SOURCE_DIR}/examples/square.points
      --tolerance 1e-12
  )
  add_test(
    NAME cli_square_k16_numerical
    COMMAND shortest_gkz
      --polygon ${CMAKE_CURRENT_SOURCE_DIR}/examples/unit_square.polygon
      --k 16
      --tolerance 1e-10
      --max-iterations 400
      --no-exact
  )
  add_test(
    NAME cli_plot_data_k2
    COMMAND shortest_gkz
      --polygon ${CMAKE_CURRENT_SOURCE_DIR}/examples/unit_square.polygon
      --k 2
      --plot-prefix ${CMAKE_CURRENT_BINARY_DIR}/test_square_k2
  )
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(
      NAME plot_results
      COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/plot_results_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_results.py
    )
    add_test(
      NAME plot_iterations
      COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/plot_iterations_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_iterations.py
    )
    add_test(
      NAME generate_wang_zhou
      COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/generate_wang_zhou_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/generate_wang_zhou.py
        $<TARGET_FILE:shortest_gkz>
    )
  endif()
endif()

# K-stability：ell_P 精确计算与相对 K-不稳定性检测（自包含子目录）。
add_subdirectory(K-stability)

option(BUILD_K_STABILITY_HIGHDIM
  "Build the independent dynamic-dimensional K-stability project" ON)
if(BUILD_K_STABILITY_HIGHDIM)
  add_subdirectory(K-stability_highdim)
endif()
