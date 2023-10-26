# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(SIMPLE_BVH_TOPLEVEL_PROJECT OFF)
else()
    set(SIMPLE_BVH_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.18.0")
if(SIMPLE_BVH_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build Simple BVH is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/SimpleBVHOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/SimpleBVHOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/SimpleBVHOptions.cmake)
endif()

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  option(SIMPLE_BVH_WITH_CCACHE "Enable ccache when building IPC Toolkit" ${SIMPLE_BVH_TOPLEVEL_PROJECT})
else()
  option(SIMPLE_BVH_WITH_CCACHE "Enable ccache when building IPC Toolkit" OFF)
endif()
if(SIMPLE_BVH_WITH_CCACHE AND CCACHE_PROGRAM)
  message(STATUS "Enabling Ccache support (${CCACHE_PROGRAM})")
  set(ccacheEnv
    CCACHE_BASEDIR=${CMAKE_BINARY_DIR}
    CCACHE_SLOPPINESS=clang_index_store,include_file_ctime,include_file_mtime,locale,pch_defines,time_macros
  )
  foreach(lang IN ITEMS C CXX)
    set(CMAKE_${lang}_COMPILER_LAUNCHER
      ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM}
    )
  endforeach()
endif()

################################################################################

project(SimpleBVH
        DESCRIPTION "A simple BVH data structure."
        LANGUAGES CXX)

option(SIMPLE_BVH_BUILD_TESTS "Build unit-tests" ${SIMPLE_BVH_TOPLEVEL_PROJECT})

# Set default minimum C++ standard
if(SIMPLE_BVH_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/simple_bvh/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")

# General CMake utils
include(simple_bvh_cpm_cache)

################################################################################
# Main library
################################################################################

add_library(simple_bvh
    src/SimpleBVH/BVH.cpp
    src/SimpleBVH/Morton.cpp
)
target_include_directories(simple_bvh PUBLIC src)
add_library(simple_bvh::simple_bvh ALIAS simple_bvh)

################################################################################
# Optional Definitions
################################################################################

################################################################################
# Required Libraries
################################################################################

include(eigen)
target_link_libraries(simple_bvh PUBLIC Eigen3::Eigen)

# Extra warnings (link last for highest priority)
include(simple_bvh_warnings)
target_link_libraries(simple_bvh PRIVATE simple_bvh::warnings)

################################################################################
# Compiler options
################################################################################

# Use C++14
target_compile_features(simple_bvh PUBLIC cxx_std_14)

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(SIMPLE_BVH_TOPLEVEL_PROJECT AND SIMPLE_BVH_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()