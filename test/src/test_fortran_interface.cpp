#include <gtest/gtest.h>
#include <qpalm_fortran.h>

TEST(TestFortranInterface, unscales_solution) {
    constexpr f_int n    = 1;
    constexpr f_int m    = 1;
    constexpr f_int h_ne = 1;
    constexpr f_int a_ne = 1;
    f_int H_ptr[]        = {1, 2};
    f_int H_row[]        = {1};
    f_float H_val[]      = {1.0};
    f_float g[]          = {0.0};
    f_int A_ptr[]        = {1, 2};
    f_int A_row[]        = {1};
    f_float A_val[]      = {100.0};
    f_float c_l[]        = {200.0};
    f_float c_u[]        = {200.0};
    f_float x[n]         = {};
    f_float y[m]         = {};
    QPALMSettings settings;
    QPALMInfo info;

    qpalm_set_default_settings(&settings);
    settings.verbose = 0;
    settings.scaling = 10;
    settings.eps_abs = 1e-9;
    settings.eps_rel = 1e-9;

    qpalm_fortran_c(n, m, h_ne, H_ptr, H_row, H_val, g, 0.0, a_ne, A_ptr, A_row,
                    A_val, c_l, c_u, settings, x, y, &info);

    EXPECT_EQ(info.status_val, QPALM_SOLVED);
    EXPECT_NEAR(x[0], 2.0, 1e-6);
    EXPECT_NEAR(y[0], -0.02, 1e-6);
    EXPECT_NEAR(info.objective, 2.0, 1e-6);
}
