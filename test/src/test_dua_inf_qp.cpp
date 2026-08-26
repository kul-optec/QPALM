#include "minunit.h"
#include <qpalm.h>

#define N 2
#define M 3
#define ANZMAX 6
#define QNZMAX 2

// Structures
static QPALMWorkspace *work; // Workspace
static QPALMSettings *settings;
static QPALMData *data;
static solver_common *c;
static solver_common common;

void dua_inf_qp_suite_setup(void) {
    settings = (QPALMSettings *)qpalm_malloc(sizeof(QPALMSettings));

    data = (QPALMData *)qpalm_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;

    data->q = (c_float *)qpalm_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = -2; 
    data->bmin = (c_float *)qpalm_calloc(M,sizeof(c_float));
    data->bmin[0] = -500; data->bmin[1] = -0.1; data->bmin[2] = -200;
    data->bmax = (c_float *)qpalm_calloc(M,sizeof(c_float));
    data->bmax[0] = 500; data->bmax[1] = 0.1; data->bmax[2] = 200;
    
    c = &common;
    solver_sparse *A = ladel_sparse_alloc(M, N, ANZMAX, UNSYMMETRIC, TRUE, FALSE);
    solver_sparse *Q = ladel_sparse_alloc(N, N, QNZMAX, UPPER, TRUE, FALSE);

    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 100.0; Ax[1] = 0.01; Ax[2] = 10.0;
    Ax[3] = 100.0; Ax[4] = 0.01; Ax[5] = 10.0;
    Ap[0] = 0; Ap[1] = 3; Ap[2] = 6;
    Ai[0] = 0; Ai[1] = 1; Ai[2] = 2; Ai[3] = 0; Ai[4] = 1; Ai[5] = 2;


    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 1e-10; Qx[1] = 1e-10; 
    Qp[0] = 0; Qp[1] = 1; Qp[2] = 2;
    Qi[0] = 0; Qi[1] = 1; 

    data->A = A;
    data->Q = Q;
}

void dua_inf_qp_suite_teardown(void) {
    qpalm_free(settings);
    // Clean setup
    data->Q = ladel_sparse_free(data->Q);
    data->A = ladel_sparse_free(data->A);
    qpalm_free(data->q);
    qpalm_free(data->bmin);
    qpalm_free(data->bmax);
    qpalm_free(data);
}


void dua_inf_qp_test_teardown(void) {
    qpalm_cleanup(work);
}

struct TestDuaInfQP : ::testing::TestWithParam<int> {
    void SetUp() override {
        dua_inf_qp_suite_setup();
        // dua_inf_qp_test_setup();
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->factorization_method = GetParam();
    settings->max_rank_update_fraction = 1.0;
    }
    void TearDown() override {
        dua_inf_qp_test_teardown();
        dua_inf_qp_suite_teardown();
    }
};

static void expect_valid_dual_infeasibility_certificate() {
    const c_float *v = work->delta_x;
    const c_float norm_v = c_max(c_absval(v[0]), c_absval(v[1]));
    const c_float qv = data->q[0] * v[0] + data->q[1] * v[1];
    const c_float norm_Av = 100 * c_absval(v[0] + v[1]);

    ASSERT_GT(norm_v, 0);
    EXPECT_LE(c_max(c_absval(data->Q->x[0] * v[0]),
                    c_absval(data->Q->x[1] * v[1])),
              settings->eps_dual_inf * norm_v);
    EXPECT_LE(qv, -settings->eps_dual_inf * norm_v);
    EXPECT_LE(norm_Av, settings->eps_dual_inf * norm_v);
}

TEST_P(TestDuaInfQP, test_dua_inf_qp) {
    settings->proximal = TRUE;
    settings->scaling = 2;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
    expect_valid_dual_infeasibility_certificate();
}
TEST_P(TestDuaInfQP, test_dua_inf_qp_unscaled) {
    settings->proximal = TRUE;
    settings->scaling = 0;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
    expect_valid_dual_infeasibility_certificate();
}
TEST_P(TestDuaInfQP, test_dua_inf_qp_noprox) {
    //This will crash actually, hence the large gamma value
    // settings->proximal = FALSE;
    settings->proximal = FALSE;
    // settings->gamma_init = 1e13;
    // settings->gamma_max = settings->gamma_init;
    settings->scaling = 2;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
    expect_valid_dual_infeasibility_certificate();
}
TEST_P(TestDuaInfQP, test_dua_inf_qp_noprox_unscaled) {
    //This will crash actually, hence the large gamma value
    // settings->proximal = FALSE;
    settings->proximal = FALSE;
    // settings->gamma_init = 1e13;
    // settings->gamma_max = settings->gamma_init;
    settings->scaling = 0;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
    expect_valid_dual_infeasibility_certificate();
}

INSTANTIATE_TEST_SUITE_P(TestDuaInfQP, TestDuaInfQP,
                         ::testing::Values(FACTORIZE_KKT_OR_SCHUR, 
                                           FACTORIZE_KKT,
                                           FACTORIZE_SCHUR));
