# src/CMakeLists.txt
#
# Build the Join-Sampling-HighDims core library (target: sjs).
#
# Expected usage from repository root:
#   add_subdirectory(src)
#
# This target is linked by CLI apps under ./apps and tests under ./tests.

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

# --------------------------
# Library target
# --------------------------

# Collect all .cpp files under src/.
file(GLOB_RECURSE SJS_LIB_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_library(sjs STATIC ${SJS_LIB_SOURCES})

target_compile_features(sjs PUBLIC cxx_std_17)

# Public headers live in ${SJS_ROOT_DIR}/include
# Some internal-only headers may live in ${SJS_ROOT_DIR}/src (rare, but allowed)
target_include_directories(sjs
  PUBLIC
    "${SJS_ROOT_DIR}/include"
    "${SJS_ROOT_DIR}/src"
)

# Reasonable warnings
if (MSVC)
  target_compile_options(sjs PRIVATE /W4 /permissive-)
else()
  target_compile_options(sjs PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Threads are used by logging and optional parallelism.
find_package(Threads REQUIRED)
target_link_libraries(sjs PUBLIC Threads::Threads)

# std::filesystem linking workaround for older GCC (<9).
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(sjs PUBLIC stdc++fs)
  endif()
endif()
