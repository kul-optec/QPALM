// THIS VERSION: 25/04/2022 AT 13:45 GMT
// Nick Gould (nick.gould@stfc.ac.uk)

#include "qpalm_fortran.h"

#define TRUE 1
#define FALSE 0

void qpalm_fortran_c( f_int n,
                      f_int m,
                      f_int h_ne,
                      f_int H_ptr[],
                      f_int H_row[],
                      f_float H_val[],
                      f_float g[n],
                      f_float f,
                      f_int a_ne,
                      f_int A_ptr[],
                      f_int A_row[],
                      f_float A_val[],
                      f_float c_l[],
                      f_float c_u[],
                      QPALMSettings settings_c,
                      f_float x[],
                      f_float y[],
                      QPALMInfo *info_c ) {

  // Structures
  QPALMWorkspace *work; // Workspace
  QPALMData *data;      // QPALMData
  
  #ifdef QPALM_FORTRAN_DEBUG_PRINT
  qpalm_print("m = %li \n",m);
  qpalm_print("n = %li \n",n);
  #endif

  // Populate QP data
  data = (QPALMData *)qpalm_malloc(sizeof(QPALMData));
  data->n = n;
  data->m = m;
  data->c = f;
  data->q = (c_float *)qpalm_malloc(n * sizeof(c_float));
  data->bmin = (c_float *)qpalm_malloc(m * sizeof(c_float));
  data->bmax = (c_float *)qpalm_malloc(m * sizeof(c_float));

  for (int i = 0; i < n; i++) {
    data->q[i] = g[i];
  }
  for (int i = 0; i < m; i++) {
    data->bmin[i] = c_l[i];
    data->bmax[i] = c_u[i];
  }

  // Populate A data
  ladel_int a_nrow;
  ladel_int a_ncol; 
  ladel_int a_nzmax; 
  ladel_int a_symmetry;

  a_nrow = m;
  a_ncol = n; 
  a_nzmax = a_ne; 
  a_symmetry = 0;

  #ifdef QPALM_FORTRAN_DEBUG_PRINT
  qpalm_print("nrow = %li \n",a_nrow);
  qpalm_print("ncol = %li \n",a_ncol);
  qpalm_print("nzmax = %li \n",a_nzmax);
  qpalm_print("\nset A\n");
  #endif

  ladel_sparse_matrix *A;

  A = ladel_sparse_alloc( a_nrow, a_ncol, a_nzmax, a_symmetry, TRUE, FALSE );

  #ifdef QPALM_FORTRAN_DEBUG_PRINT
  qpalm_print("nrow = %li \n",A->nrow);
  qpalm_print("ncol = %li \n",A->ncol);
  qpalm_print("nzmax = %li \n",A->nzmax);
  qpalm_print("symmetry = %li \n",A->symmetry);
  qpalm_print("values = %li \n",A->values);
  #endif

  c_float *Ax;
  c_int *Ai, *Ap;

  Ax = A->x;
  Ap = A->p;
  Ai = A->i;

  // N.B. convert to 0-based indices
  for (int i = 0; i < n+1; i++) {
    Ap[i] = A_ptr[i]-1;
    #ifdef QPALM_FORTRAN_DEBUG_PRINT
    printf("A column pointer %i = %i \n",i, A->p[i]);
    #endif
  }
  for (int i = 0; i < a_ne; i++) {
    Ai[i] = A_row[i]-1;
    #ifdef QPALM_FORTRAN_DEBUG_PRINT
    printf("A row index %i = %i \n",i, A->i[i]);
    #endif
  }
  for (int i = 0; i < a_ne; i++) {
    Ax[i] = A_val[i];
    #ifdef QPALM_FORTRAN_DEBUG_PRINT
    printf("A value %i = %.10f \n",i, A->x[i]);
    #endif
  }

  data->A = A;

  // Populate Q data
  ladel_int q_nrow;
  ladel_int q_ncol; 
  ladel_int q_nzmax; 
  ladel_int q_symmetry;

  q_nrow = n;
  q_ncol = n; 
  q_nzmax = h_ne; 
  q_symmetry = 1;

  #ifdef QPALM_FORTRAN_DEBUG_PRINT
  printf("nrow = %li \n",q_nrow);
  printf("ncol = %li \n",q_ncol);
  printf("nzmax = %li \n",q_nzmax);
  printf("\nset Q\n");
  #endif

  ladel_sparse_matrix *Q;

  Q = ladel_sparse_alloc( q_nrow, q_ncol, q_nzmax, q_symmetry, TRUE, FALSE );

  #ifdef QPALM_FORTRAN_DEBUG_PRINT
  printf("nrow = %li \n",Q->nrow);
  printf("ncol = %li \n",Q->ncol);
  printf("nzmax = %li \n",Q->nzmax);
  printf("symmetry = %li \n",Q->symmetry);
  printf("values = %li \n",Q->values);
  #endif

  c_float *Qx;
  c_int *Qi, *Qp;

  Qx = Q->x;
  Qp = Q->p;
  Qi = Q->i;

  // N.B. convert to 0-based indices
  for (int i = 0; i < n+1; i++) {
    Qp[i] = H_ptr[i]-1;
    #ifdef QPALM_FORTRAN_DEBUG_PRINT
    printf("Q column pointer %i = %i \n",i, Q->p[i]);
    #endif
  }
  for (int i = 0; i < h_ne; i++) {
    Qi[i] = H_row[i]-1;
    #ifdef QPALM_FORTRAN_DEBUG_PRINT
    printf("Q row index %i = %i \n",i, Q->i[i]);
    #endif
  }
  for (int i = 0; i < h_ne; i++) {
    Qx[i] = H_val[i];
    #ifdef QPALM_FORTRAN_DEBUG_PRINT
    printf("Q value %i = %.10f \n",i, Q->x[i]);
    #endif
  }

  data->Q = Q;

  // setup workspace
  work = qpalm_setup(data, &settings_c);

  // solve Problem
  qpalm_solve(work);

  // print details of solution if required
  if ( settings_c.verbose == 1 ) {
    qpalm_print("Solver status: %s\n", work->info->status);
    qpalm_print("Iter: %d\n", (int) work->info->iter);
    qpalm_print("Iter Out: %d\n", (int) work->info->iter_out);
    #ifdef QPALM_TIMING
    qpalm_print("Setup time: %f\n", work->info->setup_time);
    qpalm_print("Solve time: %f\n", work->info->solve_time);
    qpalm_print("Run time: %f\n", work->info->run_time);
    #endif
  }

  // recover solution 
  for (int i = 0; i < n; i++) {
    x[i] = work->solution->x[i];
    // printf("Solution variable %.10f \n",x[i]);
  }
  for (int i = 0; i < m; i++) {
    y[i] = work->solution->y[i];
    // printf("Multiplier variable %.10f \n",y[i]);
  }

  // return info
  *info_c = *work->info;

  // Clean workspace
  data->Q = ladel_sparse_free(data->Q);
  data->A = ladel_sparse_free(data->A);
  qpalm_free(data->q);
  qpalm_free(data->bmin);
  qpalm_free(data->bmax);

  qpalm_cleanup(work);
  qpalm_free(data);
}
