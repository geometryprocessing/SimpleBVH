################################################################################
# Tests
################################################################################

add_executable(simple_bvh_tests
    main.cpp
    tests.cpp
)

################################################################################
# Required Libraries
################################################################################

target_link_libraries(simple_bvh_tests PUBLIC simple_bvh::simple_bvh)

include(simple_bvh_warnings)
target_link_libraries(simple_bvh_tests PRIVATE simple_bvh::warnings)

include(catch2)
target_link_libraries(simple_bvh_tests PUBLIC Catch2::Catch2)

################################################################################
# Compiler options
################################################################################

################################################################################
# Register tests
################################################################################

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(simple_bvh_tests)
