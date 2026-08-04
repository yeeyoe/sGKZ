# K-stability：ell_P 精确计算与相对 K-不稳定性检测。
# 自包含子目录：只依赖 CGAL（Gmpq），不链接 gkz_core。

if(NOT TARGET CGAL::CGAL)
  find_package(CGAL REQUIRED)
endif()

add_executable(k_stability
  main.cpp
  k_stability.cpp
)
target_include_directories(k_stability PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(k_stability PRIVATE CGAL::CGAL)
target_compile_options(k_stability PRIVATE -Wall -Wextra -Wpedantic)

if(BUILD_TESTING)
  add_executable(k_stability_tests
    k_stability_tests.cpp
    k_stability.cpp
  )
  target_include_directories(k_stability_tests
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(k_stability_tests PRIVATE CGAL::CGAL)
  target_compile_options(k_stability_tests PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME k_stability_tests COMMAND k_stability_tests)

  # 单位正方形半稳定：期望退出码 2（no_counterexample_found），
  # 且输出正确的 ell_P。用 sh 包装以同时检查退出码与输出。
  add_test(
    NAME k_stability_cli_square
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability> --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon > ${CMAKE_CURRENT_BINARY_DIR}/square.out 2>&1; rc=$?; test $rc -eq 2 && grep -q 'ell_P(x,y) = 4' ${CMAKE_CURRENT_BINARY_DIR}/square.out"
  )

  add_test(
    NAME k_stability_cli_square_check_line
    COMMAND k_stability
      --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon
      --check-line "1 0 -1/2"
  )
  set_tests_properties(k_stability_cli_square_check_line PROPERTIES
    PASS_REGULAR_EXPRESSION "M_l=1/4")

  add_test(
    NAME k_stability_cli_square_svg
    COMMAND /bin/sh -c "$<TARGET_FILE:k_stability> --polygon ${CMAKE_SOURCE_DIR}/examples/unit_square.polygon --svg ${CMAKE_CURRENT_BINARY_DIR}/unit_square.svg > /dev/null 2>&1; rc=$?; test $rc -eq 2 && grep -q '</svg>' ${CMAKE_CURRENT_BINARY_DIR}/unit_square.svg"
  )
endif()
