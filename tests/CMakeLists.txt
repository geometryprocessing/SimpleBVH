cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
	main.cpp
	tests.cpp
)

ebw_download_catch2()
add_library(catch INTERFACE)
target_include_directories(catch SYSTEM INTERFACE ${EBW_EXTERNAL}/Catch2/single_include/catch2)
list(APPEND CMAKE_MODULE_PATH ${EBW_EXTERNAL}/Catch2/contrib)


add_executable(unit_tests ${test_sources})

foreach(source IN ITEMS ${test_sources})
    source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(unit_tests PUBLIC BVH_lib catch)

set(EBW_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(DATA_DIR "${CMAKE_CURRENT_LIST_DIR}/data/")
target_compile_definitions(unit_tests PUBLIC -DEBW_DATA_DIR=\"${DATA_DIR}\")



# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
