set(EBW_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(EBWDownloadExternal)

# Eigen
if(NOT TARGET Eigen3::Eigen)
    ebw_download_eigen()
    add_library(ebw_eigen INTERFACE)
    target_include_directories(ebw_eigen INTERFACE ${EBW_EXTERNAL}/eigen)
    set_property(TARGET ebw_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
    add_library(Eigen3::Eigen ALIAS ebw_eigen)
endif()



