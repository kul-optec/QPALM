cmake_minimum_required(VERSION 3.23...4.1)
set(CMAKE_CXX_SCAN_FOR_MODULES Off)
project(QPALM VERSION 1.2.6)
set(PY_VERSION_SUFFIX "")

# Options
include(CMakeDependentOption)
option(QPALM_WARNINGS_AS_ERRORS "Treat warnings as errors" Off)
option(QPALM_NONCONVEX "Compile also for nonconvex QPs" On)
option(QPALM_TIMING "Enable timing" On)
option(QPALM_PRINTING "Allow printing" On)
option(QPALM_WITH_CXX "Build the C++ interface" Off)
option(QPALM_WITH_PYTHON "Build the Python interface" Off)
option(QPALM_WITH_JULIA "Build the Julia interface" Off)
option(QPALM_WITH_MEX "Build the MATLAB interface" Off)
option(QPALM_WITH_FORTRAN "Build the Fortran interface" Off)
option(QPALM_WITH_MTX "Build the .mtx interface" Off)
option(QPALM_WITH_QPS "Build the .qps interface" Off)

cmake_dependent_option(QPALM_FORTRAN_SINGLE_PRECISION 
    "Use single-precision floating point for the Fortran interface" 
    Off QPALM_WITH_FORTRAN Off)
cmake_dependent_option(QPALM_FORTRAN_64BIT_INDICES
    "Use 64-bit integers for indices in the Fortran interface" 
    On QPALM_WITH_FORTRAN On)

include(cmake/GlobalCompileOptions.cmake)

if (QPALM_WITH_MEX AND LADEL_SINGLE_PRECISION)
    message(FATAL_ERROR "The MEX interface is not available in single-precision")
endif()

# LADEL linear algebra C library
function(add_ladel)
    set(LADEL_WITH_MEX ${QPALM_WITH_MEX})
    add_subdirectory(../LADEL/LADEL LADEL EXCLUDE_FROM_ALL)
endfunction()
find_package(LADEL QUIET)
if (NOT LADEL_FOUND)
    message(STATUS "Using included LADEL")
    add_ladel()
endif()

# QPALM C library
add_library(qpalm-headers INTERFACE)
target_include_directories(qpalm-headers
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:include>)
target_compile_definitions(qpalm-headers INTERFACE
    QPALM_USE_LADEL
    $<$<BOOL:${QPALM_NONCONVEX}>:QPALM_NONCONVEX>
    $<$<BOOL:${QPALM_TIMING}>:QPALM_TIMING>
    $<$<BOOL:${QPALM_PRINTING}>:QPALM_PRINTING>
    QPALM_VERSION_STR=\"${QPALM_VERSION}${PY_VERSION_SUFFIX}\"
    QPALM_VERSION_MAJOR=${QPALM_VERSION_MAJOR}
    QPALM_VERSION_MINOR=${QPALM_VERSION_MINOR}
    QPALM_VERSION_PATCH=${QPALM_VERSION_PATCH})
target_compile_features(qpalm-headers INTERFACE c_std_11)
if (MSVC)
    target_compile_options(qpalm-headers INTERFACE /experimental:c11atomics)
endif()
target_link_libraries(qpalm-headers INTERFACE LADEL::ladel)

add_library(qpalm-obj OBJECT
    "src/iteration.c"
    "src/lin_alg.c"
    "src/linesearch.c"
    "src/newton.c"
    "src/nonconvex.c"
    "src/qpalm.c"
    "src/scaling.c"
    "src/solver_interface.c"
    "src/termination.c"
    "src/util.c"
    "src/validate.c"
    "src/global_opts.c"
)
target_link_libraries(qpalm-obj PUBLIC qpalm-headers)
target_link_libraries(qpalm-obj PRIVATE qpalm_lax_warnings)
set_target_properties(qpalm-obj PROPERTIES C_VISIBILITY_PRESET hidden)
target_compile_definitions(qpalm-obj
    PRIVATE QPALM_EXPORTS
    INTERFACE $<$<STREQUAL:$<TARGET_PROPERTY:qpalm,TYPE>,SHARED_LIBRARY>:QPALM_IMPORTS>)
if (NOT WIN32)
    target_link_libraries(qpalm-obj PRIVATE m)
endif()

add_library(qpalm)
target_link_libraries(qpalm PUBLIC qpalm-obj)
set_target_properties(qpalm PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
add_library(${PROJECT_NAME}::qpalm ALIAS qpalm)

# C++, Python, Julia, Matlab interfaces
add_subdirectory(interfaces)

# Installation
include(cmake/Install.cmake)

# Packaging
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Pieter P")
include(CPack)
