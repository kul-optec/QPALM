cmake_minimum_required(VERSION 3.17)

add_library(qpalm_jll SHARED "qpalm_julia.c")
set_target_properties(qpalm_jll PROPERTIES C_VISIBILITY_PRESET hidden)
target_compile_definitions(qpalm_jll
    PRIVATE QPALM_EXPORTS
    INTERFACE QPALM_IMPORTS)
target_link_libraries(qpalm_jll PUBLIC qpalm-obj)
target_link_libraries(qpalm_jll PRIVATE qpalm_lax_warnings)
if (CMAKE_SYSTEM_NAME MATCHES "Linux") # Make sure LADEL is included privately
    target_link_options(qpalm_jll PRIVATE 
        "LINKER:--exclude-libs,ALL")
endif()

# Install the QPALM shared library Julia module
include(GNUInstallDirs)
install(TARGETS qpalm_jll
        COMPONENT julia_modules
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
