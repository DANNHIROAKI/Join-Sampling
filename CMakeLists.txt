# CMakeLists.txt (project root)
#
# SpatialJoinSamplingCpp: C++ experiment framework for uniform sampling from
# spatial join results (SIGMOD-style evaluation).
#
# - Library target: sjs (built from src/, headers in include/)
# - CLI apps: built from apps/
# - Tests: built from tests/
#
# This project intentionally avoids heavy dependencies (no Python, no JSON libs).
#
# Build:
#   mkdir -p build && cd build
#   cmake .. -DCMAKE_BUILD_TYPE=Release
#   cmake --build . -j
#
# Run:
#   ./sjs_run --help

cmake_minimum_required(VERSION 3.16)

project(SpatialJoinSamplingCpp LANGUAGES CXX)

# --------------------------
# Options
# --------------------------
option(SJS_BUILD_ROOT_APPS "Build command-line apps under ./apps" ON)
option(SJS_BUILD_TESTS "Build tests under ./tests (ctest)" ON)

# Most code is templated; turning off exceptions can be tricky. Keep defaults.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compilation database (nice for clangd).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --------------------------
# Build the core library
add_subdirectory(src)

if (NOT TARGET sjs)
  message(FATAL_ERROR "Target 'sjs' not found after add_subdirectory(src).")
endif()

# --------------------------
# Root apps (CLI tools)
# --------------------------
if (SJS_BUILD_ROOT_APPS)
  # Helper: consistent flags for all executables.
  function(sjs_add_app exe)
    add_executable(${exe} ${ARGN})
    target_link_libraries(${exe} PRIVATE sjs)
    target_compile_features(${exe} PRIVATE cxx_std_17)

    if (MSVC)
      target_compile_options(${exe} PRIVATE /W4 /permissive-)
    else()
      target_compile_options(${exe} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
  endfunction()

  sjs_add_app(sjs_gen_dataset  apps/sjs_gen_dataset.cpp)
  sjs_add_app(sjs_run          apps/sjs_run.cpp)
  sjs_add_app(sjs_sweep        apps/sjs_sweep.cpp)
  sjs_add_app(sjs_verify       apps/sjs_verify.cpp)
endif()

# --------------------------
# Tests
# --------------------------
if (SJS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
