# tests/CMakeLists.txt
#
# Build SJS tests.
#
# - Tests link against the "sjs" static library built in src/.
# - No third-party unit test framework is required; each test is a small
#   executable returning 0 on pass and non-zero on fail.

cmake_minimum_required(VERSION 3.16)

# Standalone mode: allow configuring tests/ directly.
# In that case, we also build the sjs library from ../src.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(sjs_tests LANGUAGES CXX)
  enable_testing()
  # Build the core library from the repository's src/.
  add_subdirectory(../src src_build)
endif()

if (NOT TARGET sjs)
  message(FATAL_ERROR "Target 'sjs' not found. Please add_subdirectory(src) before add_subdirectory(tests).")
endif()

# Helper to define a test executable with consistent flags.
function(sjs_add_test exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} PRIVATE sjs)
  target_compile_features(${exe} PRIVATE cxx_std_17)

  if (MSVC)
    target_compile_options(${exe} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${exe} PRIVATE -Wall -Wextra -Wpedantic)
  endif()

  add_test(NAME ${exe} COMMAND ${exe})
endfunction()

sjs_add_test(sjs_test_geometry         test_geometry.cpp)
sjs_add_test(sjs_test_event_sweep      test_event_sweep.cpp)
sjs_add_test(sjs_test_join_oracle      test_join_oracle.cpp)
sjs_add_test(sjs_test_sampling_quality test_sampling_quality.cpp)
sjs_add_test(sjs_test_baselines_smoke  test_baselines_smoke.cpp)
sjs_add_test(sjs_test_write_results    test_write_results.cpp)
