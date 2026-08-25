include(GNUInstallDirs)

set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/QPALM")

# Add the qpalm library to the "export-set", install the library files
install(TARGETS qpalm_fortran
    EXPORT QPALM_fortranTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        COMPONENT shlib
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" 
        COMPONENT lib)

# Install the header files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dev
    FILES_MATCHING REGEX "/.*\.[hti](pp)?$")
install(DIRECTORY "${QPALM_FORTRAN_MODULE_DIR}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dev
    FILES_MATCHING REGEX "/.*\.mod$")

# Install the export set for use with the install tree
install(EXPORT QPALM_fortranTargets 
    FILE QPALM_component_fortranTargets.cmake
    DESTINATION "${INSTALL_CMAKE_DIR}"
        COMPONENT dev
    NAMESPACE ${PROJECT_NAME}::)

# Generate the config file that includes the exports
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in"
    "${PROJECT_BINARY_DIR}/QPALM_component_fortranConfig.cmake"
    INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Install the QPALM_component_fortranConfig.cmake
install(FILES
    "${PROJECT_BINARY_DIR}/QPALM_component_fortranConfig.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}" 
        COMPONENT dev)
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/legacy-QPALM_fortranConfig.cmake"
    RENAME "QPALM_fortranConfig.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}_fortran" 
        COMPONENT dev)

# Add all targets to the build tree export set
export(EXPORT QPALM_fortranTargets
    FILE "${PROJECT_BINARY_DIR}/QPALM_component_fortranTargets.cmake"
    NAMESPACE ${PROJECT_NAME}::)
