#include "minunit.h"
#include <qpalm.h>

#define N 4
#define M 5
#define ANZMAX 4
#define QNZMAX 4

static QPALMWorkspace *work; // Workspace
static QPALMSettings *settings;
static QPALMData *data;
static solver_common *c;
static solver_common common;
static c_float solution[N] = {2.0000000e+00, -6.3801365e+01, -3.3821109e+03, -6.0483288e+00};

void basic_qp_suite_setup(void) {
    settings = (QPALMSettings *)qpalm_malloc(sizeof(QPALMSettings));

    data = (QPALMData *)qpalm_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;

    c = &common;
    data->A = ladel_sparse_alloc(M, N, ANZMAX, UNSYMMETRIC, TRUE, FALSE);
    data->Q = ladel_sparse_alloc(N, N, QNZMAX, UPPER, TRUE, FALSE);
    
    c_float *Ax = data->A->x;
    c_int *Ai = data->A->i;
    c_int *Ap = data->A->p;
    Ax[0] = -1.0000000000000000;
    Ai[0] = 3;
    Ax[1] = 0.0254311360000000;
    Ai[1] = 4;
    Ax[2] = -0.0001000000000000;
    Ai[2] = 0;
    Ax[3] = 0.3306698500000000;
    Ai[3] = 2;
    Ap[0] = 0;
    Ap[1] = 1;
    Ap[2] = 2;
    Ap[3] = 3;
    Ap[4] = 4;

    c_float *Qx = data->Q->x;
    c_int *Qi = data->Q->i;
    c_int *Qp = data->Q->p;
    Qx[0] = 1.0000000000000000;
    Qi[0] = 0;
    Qx[1] = 0.0464158880000000;
    Qi[1] = 1;
    Qx[2] = 0.0021544347000000;
    Qi[2] = 2;
    Qx[3] = 0.0001000000000000;
    Qi[3] = 3;
    Qp[0] = 0;
    Qp[1] = 1;
    Qp[2] = 2;
    Qp[3] = 3;
    Qp[4] = 4;

    data->q = (c_float *)qpalm_calloc(N,sizeof(c_float));
    data->bmin = (c_float *)qpalm_calloc(M,sizeof(c_float));
    data->bmax = (c_float *)qpalm_calloc(M,sizeof(c_float));
    data->bmin[0] = -2; 
    data->bmin[1] = -2;
    data->bmin[2] = -2;
    data->bmin[3] = -2;
    data->bmin[4] = -2;

    data->bmax[0] = 2; 
    data->bmax[1] = 2;
    data->bmax[2] = 2;
    data->bmax[3] = 2;
    data->bmax[4] = 2;

    data->q[0] = -2.0146781e+00; 
    data->q[1] = 2.9613971e+00; 
    data->q[2] = 7.2865370e+00; 
    data->q[3] = 7.8925204e+00;
}

void basic_qp_test_setup(void)
{
    settings->proximal = PROXIMAL;
    settings->scaling = SCALING;
    settings->warm_start = WARM_START;
    settings->max_iter = MAX_ITER;
    settings->inner_max_iter = INNER_MAX_ITER;
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->enable_dual_termination = ENABLE_DUAL_TERMINATION;
    settings->dual_objective_limit = DUAL_OBJECTIVE_LIMIT;
    settings->sigma_max = SIGMA_MAX;
    settings->time_limit = TIME_LIMIT;
    settings->gamma_init = 1e1;
}

void basic_qp_suite_teardown(void) {
    qpalm_free(settings);
    // Clean setup
    data->Q = ladel_sparse_free(data->Q);
    data->A = ladel_sparse_free(data->A);
    qpalm_free(data->q);
    qpalm_free(data->bmin);
    qpalm_free(data->bmax);
    qpalm_free(data);
}

void basic_qp_test_teardown(void) {
    qpalm_cleanup(work);
}

struct TestBasicQP : ::testing::TestWithParam<int> {
    void SetUp() override {
        basic_qp_suite_setup();
        basic_qp_test_setup();
        qpalm_set_default_settings(settings);
        settings->max_rank_update_fraction = 1.0;
        settings->factorization_method = GetParam();
    }
    void TearDown() override {
        basic_qp_test_teardown();
        basic_qp_suite_teardown();
    }
};

TEST_P(TestBasicQP, test_basic_qp) {
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

TEST_P(TestBasicQP, test_basic_qp_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}
TEST_P(TestBasicQP, test_basic_qp_noprox) {
    // Setup workspace
    settings->proximal = FALSE;
    settings->scaling = 2;
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}
TEST_P(TestBasicQP, test_basic_qp_noprox_unscaled) {
    // Setup workspace
    settings->proximal = FALSE;
    settings->scaling = 0;
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

TEST_P(TestBasicQP, test_basic_qp_warm_start) {
    // Setup workspace
    settings->scaling = 2;
    settings->proximal = TRUE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings);

    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.014, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

TEST_P(TestBasicQP, test_basic_qp_warm_start_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = TRUE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings);

    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.01, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

TEST_P(TestBasicQP, test_basic_qp_warm_start_noprox) {
    // Setup workspace
    settings->scaling = 2;
    settings->proximal = FALSE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings);
    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.01, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

TEST_P(TestBasicQP, test_basic_qp_warm_start_noprox_unscaled) {{
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = FALSE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings);
    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.01, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

basic_qp_test_teardown();
basic_qp_test_setup();

// TEST_P(TestBasicQP, test_basic_qp_warm_start_resolve)
{
    // Setup workspace
    work = qpalm_setup(data, settings);
    //Store initial guesses    
    c_float x[N]; 
    c_float y[M];
    prea_vec_copy(work->x, x, N);
    prea_vec_copy(work->y, y, M);
    
    // Solve Problem
    qpalm_solve(work);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    // Store solution
    c_float x_sol[N], y_sol[M];
    prea_vec_copy(work->solution->x, x_sol, N);
    prea_vec_copy(work->solution->y, y_sol, M);
    c_int iter = work->info->iter;
    
    // Warm start and resolve problem
    qpalm_warm_start(work, x, y);
    qpalm_solve(work);
    // mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_long_eq(work->info->iter, iter);
    c_float tol = 1e-12;
    for(c_int i = 0; i < N; i++) {
        mu_assert_double_eq(work->solution->x[i], x_sol[i], tol);
    }
    for(c_int i = 0; i < M; i++) {
        mu_assert_double_eq(work->solution->y[i], y_sol[i], tol);
    }

}

basic_qp_test_teardown();
basic_qp_test_setup();

// TEST_P(TestBasicQP, test_basic_qp_maxiter)
{
    // Setup workspace
    settings->max_iter = 1;

    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_MAX_ITER_REACHED);
}

basic_qp_test_teardown();
basic_qp_test_setup();

// TEST_P(TestBasicQP, test_basic_qp_inner_maxiter)
{
    // Setup workspace
    settings->eps_abs = 1e-8;
    settings->eps_rel = 1e-8;
    settings->inner_max_iter = 2;
    settings->max_iter = 10;

    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_MAX_ITER_REACHED);
}}

TEST_P(TestBasicQP, test_basic_qp_dual_objective) {
    // Setup workspace
    settings->enable_dual_termination = TRUE;

    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
    mu_assert_double_eq(work->info->objective, work->info->dual_objective, 2e-4);
}

TEST_P(TestBasicQP, test_basic_qp_dual_early_termination) {
    // Setup workspace
    settings->enable_dual_termination = TRUE;
    settings->dual_objective_limit = -1000000000.0;

    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_TERMINATED);
    mu_assert_long_eq(work->info->iter_out, 0); /* terminate on the first outer iteration */

    settings->enable_dual_termination = FALSE;
}

TEST_P(TestBasicQP, test_basic_qp_sigma_max) {
    // Setup workspace
    settings->sigma_max = 1e3;
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

TEST_P(TestBasicQP, test_basic_qp_time_limit) {
    // Setup workspace
    settings->time_limit = 0.01*1e-3;
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_TIME_LIMIT_REACHED);
}

INSTANTIATE_TEST_SUITE_P(TestBasicQP, TestBasicQP,
                         ::testing::Values(FACTORIZE_KKT_OR_SCHUR, 
                                           FACTORIZE_KKT,
                                           FACTORIZE_SCHUR));
