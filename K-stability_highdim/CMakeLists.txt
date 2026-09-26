cmake_minimum_required(VERSION 3.20)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(k_stability_highdim LANGUAGES CXX)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost 1.89 REQUIRED)

add_library(k_stability_highdim_core
  src/highdim.cpp
  src/delaunay_adapter.cpp
)
target_include_directories(k_stability_highdim_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(k_stability_highdim_core PUBLIC CGAL::CGAL Eigen3::Eigen)
target_compile_features(k_stability_highdim_core PUBLIC cxx_std_20)
target_compile_options(k_stability_highdim_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(k_stability_highdim src/main.cpp)
target_link_libraries(k_stability_highdim PRIVATE k_stability_highdim_core)
target_compile_options(k_stability_highdim PRIVATE -Wall -Wextra -Wpedantic)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(CTest)
endif()
if(BUILD_TESTING)
  add_executable(k_stability_highdim_tests tests/highdim_tests.cpp)
  target_link_libraries(k_stability_highdim_tests PRIVATE k_stability_highdim_core)
  target_compile_options(k_stability_highdim_tests PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME k_stability_highdim_tests COMMAND k_stability_highdim_tests)
  add_test(
    NAME k_stability_highdim_cli_cube5
    COMMAND k_stability_highdim
      --polytope ${CMAKE_CURRENT_SOURCE_DIR}/examples/unit_cube5.polytope
      --ell-only
  )
  set_tests_properties(k_stability_highdim_cli_cube5 PROPERTIES
    PASS_REGULAR_EXPRESSION "ell_P=10.*\\(0\\)x1")
  add_test(
    NAME k_stability_highdim_cli_simplex5
    COMMAND k_stability_highdim
      --polytope ${CMAKE_CURRENT_SOURCE_DIR}/examples/simplex5.polytope
      --ell-only
  )
  set_tests_properties(k_stability_highdim_cli_simplex5 PROPERTIES
    PASS_REGULAR_EXPRESSION "ell_P=30.*\\(0\\)x1")
endif()
