cmake_minimum_required(VERSION 3.20)
project(chowk LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

add_library(chowk_core
  src/chowk.cpp
  src/highdim.cpp
  src/delaunay_adapter.cpp
  src/sgkz.cpp
)
target_include_directories(chowk_core PUBLIC include src)
target_link_libraries(chowk_core PUBLIC CGAL::CGAL Eigen3::Eigen Boost::boost)
target_compile_options(chowk_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(chowk src/main.cpp)
target_link_libraries(chowk PRIVATE chowk_core)

include(CTest)
if(BUILD_TESTING)
  add_executable(chowk_tests tests/chowk_tests.cpp)
  target_link_libraries(chowk_tests PRIVATE chowk_core)
  add_test(NAME chowk_tests COMMAND chowk_tests)
endif()
