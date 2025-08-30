cmake_minimum_required(VERSION 3.22...4.1)
project(QPALM-mex)

# Find the Matlab MEX compiler:
find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER)

# Create a target for compiling a MEX function:
matlab_add_mex(
    NAME qpalm_mex                       # Name for use in Matlab
    SRC src/qpalm_mex.c                  # Source files to compile
    LINK_TO QPALM::qpalm LADEL::mex_util # Dependencies
)

# Install the final MEX file
install(FILES $<TARGET_FILE:qpalm_mex>
        qpalm.m
        EXCLUDE_FROM_ALL
        COMPONENT mex_interface
        DESTINATION +qpalm)
