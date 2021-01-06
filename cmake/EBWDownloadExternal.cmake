################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(EBW_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(EBW_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(ebw_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${EBW_EXTERNAL}/${name}
        DOWNLOAD_DIR ${EBW_EXTERNAL}/.cache/${name}
        QUIET
        ${EBW_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################


## Catch2 BSL 1.0 optional
function(ebw_download_catch2)
    ebw_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.4.2
    )
endfunction()


set(EBW_EIGEN_VERSION 3.3.7 CACHE STRING "Default version of Eigen used.")
function(ebw_download_eigen)
ebw_download_project(eigen
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG        ${LIBIGL_EIGEN_VERSION}
		${EBW_EIGEN_VERSION}
	)
endfunction()

