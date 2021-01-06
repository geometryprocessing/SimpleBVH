################################################################################
cmake_minimum_required(VERSION 3.13)
project(BVH_lib)
################################################################################

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(EBW_TOPLEVEL_PROJECT ON)
else()
    set(EBW_TOPLEVEL_PROJECT OFF)
endif()


set(EBW_EXTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


# Use C++11/14
include(CXXFeatures)

include(EBWDependencies)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################


################################################################################
# Main library
################################################################################

add_library(BVH_lib src/BVH.cpp src/Morton.cpp)
target_include_directories(BVH_lib PUBLIC src)

# Setup dependencies
target_compile_features(BVH_lib PUBLIC ${CXX14_FEATURES})
target_link_libraries(BVH_lib PUBLIC Eigen3::Eigen)

################################################################################
# Subdirectories
################################################################################

# Other sources
# add_subdirectory(src)

# Compile tests only if this is a top-level project
if(EBW_TOPLEVEL_PROJECT)
    # Unit tests
    include(CTest)
    # Catch2
    enable_testing()
    add_subdirectory(tests)
endif()
