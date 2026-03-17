# src/CMakeLists.txt
#
# Build the Spatial Join Sampling (SJS) C++ experiment framework.
#
# This CMakeLists is designed to be added via:
#   add_subdirectory(src)
# from the project root, but it can also serve as a standalone build entry
# if you run CMake from the src/ directory.

cmake_minimum_required(VERSION 3.16)

# If this directory is the top-level (standalone build), declare the project.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(sjs LANGUAGES CXX)
  # In a standalone build, the repository root is the parent of src/.
  set(SJS_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
else()
  # When included via add_subdirectory(src), PROJECT_SOURCE_DIR refers to repo root.
  set(SJS_ROOT_DIR "${PROJECT_SOURCE_DIR}")
endif()

## Note: CLI tools are built from the repository root (apps/).
## We intentionally do NOT build legacy src/apps entry points.

# --------------------------
# Library target
# --------------------------

# Collect all non-app .cpp files under src/.
# (Most implementations are header-only; .cpp files mainly provide factories and
# explicit Dim=2 instantiations.)
file(GLOB_RECURSE SJS_LIB_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/baselines/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/join/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/sampling/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/geometry/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/data/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp"
)

# Adaptive framework source is removed; keep a defensive filter for stale trees.
list(FILTER SJS_LIB_SOURCES EXCLUDE REGEX ".*/baselines/ours/adaptive\\.cpp$")

add_library(sjs STATIC ${SJS_LIB_SOURCES})

target_compile_features(sjs PUBLIC cxx_std_17)

# Headers live in:
#   - ${SJS_ROOT_DIR}/include   (public headers)
#   - ${SJS_ROOT_DIR}/src       (internal factory headers under src/)
#
# Note: baseline_factory_2d.h is intentionally under src/.
target_include_directories(sjs
  PUBLIC
    "${SJS_ROOT_DIR}"
    "${SJS_ROOT_DIR}/include"
    "${SJS_ROOT_DIR}/src"
)

# Reasonable warnings (tune as needed).
if (MSVC)
  target_compile_options(sjs PRIVATE /W4 /permissive-)
else()
  target_compile_options(sjs PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Threads are used by logging and optional future parallel baselines.
find_package(Threads REQUIRED)
target_link_libraries(sjs PUBLIC Threads::Threads)

# std::filesystem linking workaround for older GCC (<9).
# (Modern GCC/Clang/libc++ do not need this.)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(sjs PUBLIC stdc++fs)
  endif()
endif()
