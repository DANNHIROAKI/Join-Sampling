# CMakeLists.txt (project root)
#
# Join-Sampling-HighDims: C++ experiment framework for uniform sampling
# from spatial join results, upgraded to support runtime dimension dispatch.
#
# Targets:
#   - sjs (static library from src/)
#   - Apps (from apps/):
#       sjs_run, sjs_gen_dataset, sjs_sweep, sjs_verify
#   - Tests (from tests/):
#       sjs_test_*
#
# Build:
#   mkdir -p build/release && cd build/release
#   cmake ../.. -DCMAKE_BUILD_TYPE=Release
#   cmake --build . -j
#
# Run:
#   ./sjs_run --help

cmake_minimum_required(VERSION 3.16)

project(JoinSamplingHighDims LANGUAGES CXX)

# --------------------------
# Options
# --------------------------
option(SJS_BUILD_APPS  "Build command-line apps under ./apps" ON)
option(SJS_BUILD_TESTS "Build tests under ./tests (ctest)" ON)

# Maximum dimension to compile into the library.
#
# Higher values increase compile time due to template instantiation.
# Default: 12 (covers common HighDims experiments).
set(SJS_MAX_DIM "12" CACHE STRING "Maximum supported compile-time dimension (>=2).")

# --------------------------
# C++ standard
# --------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compilation database (useful for clangd)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optional: keep runtime outputs in build root
# so scripts can find them easily.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# --------------------------
# Build the core library
# --------------------------
add_subdirectory(src)

if (NOT TARGET sjs)
  message(FATAL_ERROR "Target 'sjs' not found after add_subdirectory(src).")
endif()

# Propagate the max-dim setting to all compilation units.
# This controls SJS_DIMS(X) via sjs/dispatch/dim_dispatch.h.
target_compile_definitions(sjs PUBLIC SJS_MAX_DIM=${SJS_MAX_DIM})

# Helper: consistent flags for executables.
function(sjs_apply_exe_defaults exe_target)
  target_link_libraries(${exe_target} PRIVATE sjs)
  target_compile_features(${exe_target} PRIVATE cxx_std_17)

  if (MSVC)
    target_compile_options(${exe_target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${exe_target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

# --------------------------
# Apps
# --------------------------
if (SJS_BUILD_APPS)
  add_executable(sjs_run          apps/sjs_run.cpp)
  sjs_apply_exe_defaults(sjs_run)

  add_executable(sjs_gen_dataset  apps/sjs_gen_dataset.cpp)
  sjs_apply_exe_defaults(sjs_gen_dataset)

  add_executable(sjs_sweep        apps/sjs_sweep.cpp)
  sjs_apply_exe_defaults(sjs_sweep)

  add_executable(sjs_verify       apps/sjs_verify.cpp)
  sjs_apply_exe_defaults(sjs_verify)
endif()

# --------------------------
# Tests
# --------------------------
if (SJS_BUILD_TESTS)
  enable_testing()

  function(sjs_add_test exe_name src_file)
    add_executable(${exe_name} ${src_file})
    sjs_apply_exe_defaults(${exe_name})
    add_test(NAME ${exe_name} COMMAND ${exe_name})
  endfunction()

  sjs_add_test(sjs_test_geometry_nd          tests/test_geometry_nd.cpp)
  sjs_add_test(sjs_test_event_sweep_nd       tests/test_event_sweep_nd.cpp)
  sjs_add_test(sjs_test_rank_space           tests/test_rank_space.cpp)
  sjs_add_test(sjs_test_segtree_1d           tests/test_segtree_1d.cpp)
  sjs_add_test(sjs_test_mode_index_nd        tests/test_mode_index_nd.cpp)
  sjs_add_test(sjs_test_mode_family_nd       tests/test_mode_family_nd.cpp)
  sjs_add_test(sjs_test_sampling_quality_nd  tests/test_sampling_quality_nd.cpp)
  sjs_add_test(sjs_test_baselines_smoke_nd   tests/test_baselines_smoke_nd.cpp)
endif()
