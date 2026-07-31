# CMake generated Testfile for 
# Source directory: /Users/yaoy/Documents/sGKZ
# Build directory: /Users/yaoy/Documents/sGKZ/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[gkz_tests]=] "/Users/yaoy/Documents/sGKZ/build/gkz_tests")
set_tests_properties([=[gkz_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;25;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[oracle_kernel_consistency]=] "/Users/yaoy/Documents/sGKZ/build/oracle_kernel_consistency_test")
set_tests_properties([=[oracle_kernel_consistency]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;29;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[cli_square]=] "/Users/yaoy/Documents/sGKZ/build/shortest_gkz" "--points" "/Users/yaoy/Documents/sGKZ/examples/square.points" "--tolerance" "1e-12")
set_tests_properties([=[cli_square]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;31;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[cli_square_k16_numerical]=] "/Users/yaoy/Documents/sGKZ/build/shortest_gkz" "--polygon" "/Users/yaoy/Documents/sGKZ/examples/unit_square.polygon" "--k" "16" "--tolerance" "1e-10" "--max-iterations" "400" "--no-exact")
set_tests_properties([=[cli_square_k16_numerical]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;37;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[cli_plot_data_k2]=] "/Users/yaoy/Documents/sGKZ/build/shortest_gkz" "--polygon" "/Users/yaoy/Documents/sGKZ/examples/unit_square.polygon" "--k" "2" "--plot-prefix" "/Users/yaoy/Documents/sGKZ/build/test_square_k2")
set_tests_properties([=[cli_plot_data_k2]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;46;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[plot_results]=] "/opt/homebrew/Frameworks/Python.framework/Versions/3.14/bin/python3.14" "/Users/yaoy/Documents/sGKZ/tests/plot_results_test.py" "/Users/yaoy/Documents/sGKZ/plot_results.py")
set_tests_properties([=[plot_results]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;55;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[plot_iterations]=] "/opt/homebrew/Frameworks/Python.framework/Versions/3.14/bin/python3.14" "/Users/yaoy/Documents/sGKZ/tests/plot_iterations_test.py" "/Users/yaoy/Documents/sGKZ/plot_iterations.py")
set_tests_properties([=[plot_iterations]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;61;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
add_test([=[generate_wang_zhou]=] "/opt/homebrew/Frameworks/Python.framework/Versions/3.14/bin/python3.14" "/Users/yaoy/Documents/sGKZ/tests/generate_wang_zhou_test.py" "/Users/yaoy/Documents/sGKZ/generate_wang_zhou.py" "/Users/yaoy/Documents/sGKZ/build/shortest_gkz")
set_tests_properties([=[generate_wang_zhou]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/yaoy/Documents/sGKZ/CMakeLists.txt;67;add_test;/Users/yaoy/Documents/sGKZ/CMakeLists.txt;0;")
