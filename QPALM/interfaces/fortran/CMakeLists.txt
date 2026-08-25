enable_language(Fortran)

add_library(qpalm_fortran
    "qpalm_fiface.F90"
    "qpalm_fortran.c"
    "include/qpalm_fortran.h"
)
target_link_libraries(qpalm_fortran PUBLIC QPALM::qpalm)
set(QPALM_FORTRAN_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
target_include_directories(qpalm_fortran
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
           $<INSTALL_INTERFACE:include>)
target_include_directories(qpalm_fortran
    PUBLIC $<BUILD_INTERFACE:${QPALM_FORTRAN_MODULE_DIR}>
           $<INSTALL_INTERFACE:include>)
set_target_properties(qpalm_fortran PROPERTIES
    Fortran_PREPROCESS On
    Fortran_MODULE_DIRECTORY ${QPALM_FORTRAN_MODULE_DIR})
target_compile_definitions(qpalm_fortran PUBLIC
    $<$<BOOL:${QPALM_FORTRAN_SINGLE_PRECISION}>:QPALM_FORTRAN_SINGLE_PRECISION>
    $<$<BOOL:${QPALM_FORTRAN_64BIT_INDICES}>:QPALM_FORTRAN_64BIT_INDICES>)
set_target_properties(qpalm_fortran PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
add_library(QPALM::qpalm_fortran ALIAS qpalm_fortran)

include(cmake/Install.cmake)