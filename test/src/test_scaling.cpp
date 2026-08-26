#include <gtest/gtest.h>
#include <qpalm.h>

class TestScaling : public ::testing::Test {
  protected:
    void SetUp() override {
        qpalm_set_default_settings(&settings);
        settings.verbose                 = 0;
        settings.scaling                 = 10;
        settings.eps_abs                 = 1e-9;
        settings.eps_rel                 = 1e-9;
        settings.enable_dual_termination = 1;

        data.n    = 1;
        data.m    = 1;
        data.c    = 0;
        data.q    = q;
        data.bmin = bmin;
        data.bmax = bmax;
        data.Q    = ladel_sparse_alloc(1, 1, 1, UPPER, TRUE, FALSE);
        data.A    = ladel_sparse_alloc(1, 1, 1, UNSYMMETRIC, TRUE, FALSE);

        data.Q->p[0] = 0;
        data.Q->p[1] = 1;
        data.Q->i[0] = 0;
        data.Q->x[0] = 1;
        data.A->p[0] = 0;
        data.A->p[1] = 1;
        data.A->i[0] = 0;
        data.A->x[0] = 100;
    }

    void TearDown() override {
        if (work)
            qpalm_cleanup(work);
        data.Q = ladel_sparse_free(data.Q);
        data.A = ladel_sparse_free(data.A);
    }

    QPALMSettings settings;
    QPALMData data;
    QPALMWorkspace *work = nullptr;
    c_float q[1]         = {1000};
    c_float bmin[1]      = {200};
    c_float bmax[1]      = {200};
};

TEST_F(TestScaling, warm_start_and_solution_round_trip) {
    c_float x[] = {2};
    c_float y[] = {-10.02};

    work = qpalm_setup(&data, &settings);
    ASSERT_NE(work, nullptr);
    qpalm_warm_start(work, nullptr, nullptr);
    qpalm_warm_start(work, x, y);
    qpalm_solve(work);

    ASSERT_EQ(work->info->status_val, QPALM_SOLVED);
    EXPECT_EQ(work->info->iter, 0);
    EXPECT_NE(work->scaling->D[0], 1);
    EXPECT_NE(work->scaling->E[0], 1);
    EXPECT_NE(work->scaling->c, 1);
    EXPECT_NEAR(work->yh[0], work->scaling->c * work->scaling->Einv[0] * y[0],
                1e-12);
    EXPECT_NEAR(work->solution->x[0], x[0], 1e-8);
    EXPECT_NEAR(work->solution->y[0], y[0], 1e-8);
    EXPECT_NEAR(work->info->objective, 2002, 1e-8);
    EXPECT_NEAR(work->info->dual_objective, 2002, 1e-8);

    EXPECT_DOUBLE_EQ(work->data->Q->x[0], data.Q->x[0]);
    EXPECT_DOUBLE_EQ(work->data->A->x[0], data.A->x[0]);
    EXPECT_DOUBLE_EQ(work->data->q[0], data.q[0]);
    EXPECT_DOUBLE_EQ(work->data->bmin[0], data.bmin[0]);
    EXPECT_DOUBLE_EQ(work->data->bmax[0], data.bmax[0]);

    qpalm_warm_start(work, work->solution->x, work->solution->y);
    qpalm_solve(work);

    ASSERT_EQ(work->info->status_val, QPALM_SOLVED);
    EXPECT_EQ(work->info->iter, 0);
    EXPECT_NEAR(work->solution->x[0], x[0], 1e-8);
    EXPECT_NEAR(work->solution->y[0], y[0], 1e-8);
    EXPECT_NEAR(work->info->dual_objective, 2002, 1e-8);
    EXPECT_DOUBLE_EQ(work->data->Q->x[0], data.Q->x[0]);
    EXPECT_DOUBLE_EQ(work->data->A->x[0], data.A->x[0]);
    EXPECT_DOUBLE_EQ(work->data->q[0], data.q[0]);
    EXPECT_DOUBLE_EQ(work->data->bmin[0], data.bmin[0]);
    EXPECT_DOUBLE_EQ(work->data->bmax[0], data.bmax[0]);
}

TEST_F(TestScaling, cancellation_returns_unscaled_warm_start) {
    c_float x[] = {2};
    c_float y[] = {-10.02};

    work = qpalm_setup(&data, &settings);
    ASSERT_NE(work, nullptr);
    qpalm_warm_start(work, x, y);
    qpalm_cancel(work);
    qpalm_solve(work);

    ASSERT_EQ(work->info->status_val, QPALM_USER_CANCELLATION);
    EXPECT_NEAR(work->solution->x[0], x[0], 1e-8);
    EXPECT_NEAR(work->solution->y[0], y[0], 1e-8);
    EXPECT_DOUBLE_EQ(work->data->Q->x[0], data.Q->x[0]);
    EXPECT_DOUBLE_EQ(work->data->A->x[0], data.A->x[0]);
    EXPECT_DOUBLE_EQ(work->data->q[0], data.q[0]);
    EXPECT_DOUBLE_EQ(work->data->bmin[0], data.bmin[0]);
    EXPECT_DOUBLE_EQ(work->data->bmax[0], data.bmax[0]);
}
