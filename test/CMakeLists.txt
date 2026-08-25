include(GoogleTest)
find_package(GTest REQUIRED CONFIG)

add_executable(tests
    "src/test_ls_qp.cpp"
    "src/test_lin_alg.cpp"
    "src/test_update_after_setup.cpp"
    "src/test_medium_qp.cpp"
    "src/test_update.cpp"
    "src/test_casadi_general_unconstrained.cpp"
    "src/test_dua_inf_qp.cpp"
    "src/test_validate.cpp"
    "src/test_error_handling.cpp"
    "src/test_nonconvex_qp.cpp"
    "src/test_prim_inf_qp.cpp"
    "src/test_solver_interface.cpp"
    "src/test_casadi_general_convex_sparse.cpp"
    "src/test_basic_qp.cpp"
    "src/test_degen_hess.cpp"
)
if (QPALM_WITH_FORTRAN)
    target_sources(tests PRIVATE "src/test_fortran_interface.cpp")
    target_link_libraries(tests PRIVATE QPALM::qpalm_fortran)
endif()
target_include_directories(tests
    PRIVATE include
            include/Minunit
)

target_link_libraries(tests PRIVATE qpalm-obj)
target_link_libraries(tests PRIVATE GTest::gtest_main)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.21) 
    # Copy DLLs for test executable
    add_custom_command(TARGET tests PRE_LINK COMMAND
        COMMAND ${CMAKE_COMMAND} -E 
            $<IF:$<BOOL:$<TARGET_RUNTIME_DLLS:tests>>,copy,true>
            $<TARGET_RUNTIME_DLLS:tests> $<TARGET_FILE_DIR:tests>
        COMMAND_EXPAND_LISTS)
endif()

option(QPALM_FORCE_TEST_DISCOVERY Off)
if (NOT CMAKE_CROSSCOMPILING OR QPALM_FORCE_TEST_DISCOVERY)
    gtest_discover_tests(tests DISCOVERY_TIMEOUT 60)
endif()
add_executable(QPALM::tests ALIAS tests)
