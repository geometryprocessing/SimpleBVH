################################################################################
# Tests
################################################################################

add_executable(simple_bvh_tests tests.cpp)

################################################################################
# Required Libraries
################################################################################

target_link_libraries(simple_bvh_tests PUBLIC simple_bvh::simple_bvh)

include(catch2)
target_link_libraries(simple_bvh_tests PUBLIC Catch2::Catch2WithMain)

# Extra warnings (link last for highest priority)
include(simple_bvh_warnings)
target_link_libraries(simple_bvh_tests PRIVATE simple_bvh::warnings)

################################################################################
# Compiler options
################################################################################

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(simple_bvh_tests)
