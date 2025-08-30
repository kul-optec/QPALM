cmake_minimum_required(VERSION 3.20)
set(PY_FULL_VERSION ${PROJECT_VERSION}${PY_VERSION_SUFFIX})
set(PY_BUILD_CMAKE_MODULE_NAME "qpalm" CACHE STRING "Python module name")
set(PY_BUILD_CMAKE_PACKAGE_NAME "qpalm" CACHE STRING "Python package name")

# Make sure that the Python and CMake versions match
if (DEFINED VERIFY_VERSION)
    if (NOT "${VERIFY_VERSION}" MATCHES "^${PY_FULL_VERSION}$")
        message(FATAL_ERROR "Version number does not match "
                             "(${VERIFY_VERSION} - ${PY_FULL_VERSION}).")
    endif()
endif()

# Python interface requires the C++ interface
if (NOT QPALM_WITH_CXX)
    message(FATAL_ERROR "QPALM Python interface requires the C++ interface."
                        "Enable it using -DQPALM_WITH_CXX=On.")
endif()

# Find the Python development files
if (CMAKE_CROSSCOMPILING AND NOT DEFINED CMAKE_CROSSCOMPILING_EMULATOR)
    find_package(Python3 REQUIRED COMPONENTS Development.Module)
else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
endif()

# Find the pybind11 library
set(PYBIND11_USE_CROSSCOMPILING On) # https://github.com/pybind/pybind11/pull/5083
find_package(pybind11 REQUIRED)

# Compile the example Python module
pybind11_add_module(_qpalm MODULE "qpalm.py.cpp" "demangled-typename.cpp")
target_compile_features(_qpalm PRIVATE cxx_std_17)
target_link_libraries(_qpalm PRIVATE qpalm_warnings)
target_link_libraries(_qpalm PRIVATE pybind11::pybind11 qpalm qpalm_cxx)
target_compile_definitions(_qpalm PRIVATE
    MODULE_NAME=$<TARGET_FILE_BASE_NAME:_qpalm>
    VERSION_INFO="${PY_FULL_VERSION}")
set_target_properties(_qpalm PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_options(_qpalm PRIVATE "LINKER:--exclude-libs,ALL")
endif()
set(QPALM_PYTHON_DEBUG_CONFIG "Debug" CACHE STRING "")
set_property(CACHE QPALM_PYTHON_DEBUG_CONFIG PROPERTY STRINGS
    ${CMAKE_CONFIGURATION_TYPES})
if (QPALM_PYTHON_DEBUG_CONFIG)
    string(TOUPPER ${QPALM_PYTHON_DEBUG_CONFIG} QPALM_PYTHON_DEBUG_CONFIG_UP)
    set_target_properties(_qpalm PROPERTIES
        ${QPALM_PYTHON_DEBUG_CONFIG_UP}_POSTFIX "_d"
        PDB_NAME_${QPALM_PYTHON_DEBUG_CONFIG_UP} "_qpalm_d")
endif()

if (DEFINED PY_BUILD_CMAKE_VERSION)
    set(QPALM_MODULE_INSTALL_DIR "")
else()
    set(QPALM_MODULE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
endif()

# Install the Python source code
install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/python/${PY_BUILD_CMAKE_MODULE_NAME}
        COMPONENT python_source
        DESTINATION "${QPALM_MODULE_INSTALL_DIR}"
        PATTERN "__pycache__" EXCLUDE)

# Install the Python module
install(TARGETS _qpalm
        COMPONENT python_modules
        DESTINATION ./${QPALM_MODULE_INSTALL_DIR}/${PY_BUILD_CMAKE_MODULE_NAME})

# Strip and install debug information
set(DEBUG_FILE "$<TARGET_FILE_NAME:_qpalm>.debug")
if (CMAKE_STRIP AND CMAKE_OBJCOPY)
    add_custom_command(TARGET _qpalm POST_BUILD
        COMMAND "${CMAKE_STRIP}" "--only-keep-debug" "$<TARGET_FILE:_qpalm>" "-o" "${DEBUG_FILE}"
        COMMAND "${CMAKE_STRIP}" "--strip-debug" "$<TARGET_FILE:_qpalm>"
        COMMAND "${CMAKE_OBJCOPY}" "--add-gnu-debuglink=${DEBUG_FILE}" "$<TARGET_FILE:_qpalm>"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${DEBUG_FILE}"
        EXCLUDE_FROM_ALL
        CONFIGURATIONS Debug RelWithDebInfo
        COMPONENT python_modules_debug
        DESTINATION ./${QPALM_MODULE_INSTALL_DIR}/${PY_BUILD_CMAKE_MODULE_NAME})
endif()

# Generate stubs for the Python module
set(QPALM_WITH_PY_STUBS_DEFAULT On)
if (CMAKE_CROSSCOMPILING)
    set(QPALM_WITH_PY_STUBS_DEFAULT Off)
endif()
option(QPALM_WITH_PY_STUBS
    "Generate Python stub files (.pyi) for the Python module."
    ${QPALM_WITH_PY_STUBS_DEFAULT})
if (QPALM_WITH_PY_STUBS)
    include(cmake/Pybind11Stubgen.cmake)
    pybind11_stubgen(_qpalm
        COMPONENT python_stubs
        PACKAGE ${PY_BUILD_CMAKE_MODULE_NAME}
        PACKAGE_ROOT "${QPALM_MODULE_INSTALL_DIR}")
endif()
