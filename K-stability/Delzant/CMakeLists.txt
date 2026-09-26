add_executable(unstable_Delzant
  main.cpp
  unstable_delzant.cpp
)
target_link_libraries(unstable_Delzant PRIVATE k_stability_core SQLite3::SQLite3)
target_include_directories(unstable_Delzant PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(unstable_Delzant PRIVATE -Wall -Wextra -Wpedantic)

if(BUILD_TESTING)
  add_executable(unstable_delzant_tests
    tests.cpp
    unstable_delzant.cpp
  )
  target_link_libraries(unstable_delzant_tests PRIVATE k_stability_core SQLite3::SQLite3)
  target_include_directories(unstable_delzant_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_options(unstable_delzant_tests PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME unstable_delzant_tests COMMAND unstable_delzant_tests)
  add_test(
    NAME unstable_delzant_requires_max_diameter
    COMMAND /bin/sh -c "$<TARGET_FILE:unstable_Delzant> --d 3 --time-limit 1 >/dev/null 2>&1; test $? -eq 1"
  )
  add_test(
    NAME unstable_delzant_rejects_invalid_max_diameter
    COMMAND /bin/sh -c "$<TARGET_FILE:unstable_Delzant> --d 3 --time-limit 1 --max-diameter 0 >/dev/null 2>&1; test $? -eq 1"
  )
endif()
