# K-stability：ell_P 精确计算与相对 K-不稳定性检测。
# 自包含子目录：只依赖 CGAL（Gmpq），不链接 gkz_core。

if(NOT TARGET CGAL::CGAL)
  find_package(CGAL REQUIRED)
endif()

find_package(SQLite3 REQUIRED)

add_library(k_stability_core
  k_stability.cpp
  decompose.cpp
)
target_include_directories(k_stability_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(k_stability_core PUBLIC CGAL::CGAL)
target_compile_options(k_stability_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(k_stability
  main.cpp
)
target_include_directories(k_stability PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(k_stability PRIVATE k_stability_core)

add_executable(k_stability_decompose decompose_main.cpp)
target_include_directories(k_stability_decompose PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(k_stability_decompose PRIVATE k_stability_core)
target_compile_options(k_stability_decompose PRIVATE -Wall -Wextra -Wpedantic)

add_library(k_stability_search_core
  search.cpp
)
target_include_directories(k_stability_search_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(k_stability_search_core PUBLIC k_stability_core SQLite3::SQLite3)
target_compile_options(k_stability_search_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(k_stability_search search_main.cpp)
target_link_libraries(k_stability_search PRIVATE k_stability_search_core)
target_compile_options(k_stability_search PRIVATE -Wall -Wextra -Wpedantic)

if(BUILD_TESTING)
  add_executable(k_stability_tests
    k_stability_tests.cpp
  )
  target_include_directories(k_stability_tests
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(k_stability_tests PRIVATE k_stability_core)
  add_test(NAME k_stability_tests COMMAND k_stability_tests)

  # 单位正方形半稳定：期望退出码 2（no_counterexample_found），
  # 且输出正确的 ell_P。用 sh 包装以同时检查退出码与输出。
  add_test(
    NAME k_stability_cli_square
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability> --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon > ${CMAKE_CURRENT_BINARY_DIR}/square.out 2>&1; rc=$?; test $rc -eq 2 && grep -q 'ell_P(x,y) = 4' ${CMAKE_CURRENT_BINARY_DIR}/square.out"
  )

  add_executable(k_stability_search_tests search_tests.cpp)
  target_link_libraries(k_stability_search_tests PRIVATE k_stability_search_core)
  target_compile_options(k_stability_search_tests PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME k_stability_search_tests COMMAND k_stability_search_tests)

  add_test(
    NAME k_stability_cli_square_check_line
    COMMAND k_stability
      --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon
      --check-line "1 0 -1/2"
  )
  set_tests_properties(k_stability_cli_square_check_line PROPERTIES
    PASS_REGULAR_EXPRESSION "M_l=1/4")

  add_test(
    NAME k_stability_cli_square_bottom_null_measure
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability> --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square_bottom_null_measure.polygon --check-line '1 0 -1/2' > ${CMAKE_CURRENT_BINARY_DIR}/square_bottom_null.out 2>&1 && grep -q 'boundary_length_dsigma=3' ${CMAKE_CURRENT_BINARY_DIR}/square_bottom_null.out && grep -q 'null_measure_edge_count=1' ${CMAKE_CURRENT_BINARY_DIR}/square_bottom_null.out && grep -q 'ell_P(x,y) = 0 + (0) x + (6) y' ${CMAKE_CURRENT_BINARY_DIR}/square_bottom_null.out && grep -q 'M_l=1/4' ${CMAKE_CURRENT_BINARY_DIR}/square_bottom_null.out"
  )

  add_test(
    NAME k_stability_decompose_parent_not_certified
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability_decompose> --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon --theta-steps 24 --t-steps 24 > /dev/null 2>&1; test $? -eq 2"
  )

  add_test(
    NAME k_stability_decompose_d5_a74
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability_decompose> --polygon ${CMAKE_SOURCE_DIR}/examples/d5_a74 --theta-steps 180 --t-steps 128 --root-starts 4 --root-iterations 30 --root-tolerance 1e-8 > ${CMAKE_CURRENT_BINARY_DIR}/d5_a74_decompose.out 2>&1 && grep -q 'parent_certified=true' ${CMAKE_CURRENT_BINARY_DIR}/d5_a74_decompose.out && grep -q 'decomposition_candidates=1' ${CMAKE_CURRENT_BINARY_DIR}/d5_a74_decompose.out && grep -q 'candidate_concave=true' ${CMAKE_CURRENT_BINARY_DIR}/d5_a74_decompose.out"
  )

  add_test(
    NAME k_stability_cli_square_svg
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability> --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon --svg ${CMAKE_CURRENT_BINARY_DIR}/unit_square.svg > /dev/null 2>&1; rc=$?; test $rc -eq 2 && grep -q '</svg>' ${CMAKE_CURRENT_BINARY_DIR}/unit_square.svg"
  )
endif()
