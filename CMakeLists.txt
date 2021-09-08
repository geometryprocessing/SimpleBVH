# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(SIMPLE_BVH_TOPLEVEL_PROJECT OFF)
else()
    set(SIMPLE_BVH_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
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

################################################################################

project(SimpleBVH
        DESCRIPTION "A simple BVH data structure."
        LANGUAGES CXX)

option(SIMPLE_BVH_WITH_UNIT_TESTS "Build unit-tests" ${SIMPLE_BVH_TOPLEVEL_PROJECT})

# Set default minimum C++ standard
if(RIGID_IPC_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/simple_bvh/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")

################################################################################
# Main library
################################################################################

add_library(simple_bvh
    src/BVH.cpp
    src/Morton.cpp
)
target_include_directories(simple_bvh PUBLIC src)
add_library(simple_bvh::simple_bvh ALIAS simple_bvh)

################################################################################
# Optional Definitions
################################################################################

################################################################################
# Required Libraries
################################################################################

# Extra warnings
include(simple_bvh_warnings)
target_link_libraries(simple_bvh PRIVATE simple_bvh::warnings)

include(eigen)
target_link_libraries(simple_bvh PUBLIC Eigen3::Eigen)

################################################################################
# Compiler options
################################################################################

# Use C++14
target_compile_features(simple_bvh PUBLIC cxx_std_14)

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(SIMPLE_BVH_WITH_UNIT_TESTS)
    include(CTest)
    enable_testing()

    # Include Catch2 and provide function `catch_discover_tests` to register tests.
    include(catch2)
    FetchContent_GetProperties(catch2)
    include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")

    add_subdirectory(tests)
endif()
