/**
 * @file scaling.c
 * @author Ben Hermans
 * @brief Problem data scaling during setup.
 * @details This file includes the routine that is called during setup to scale the problem data, 
 * and initial guesses if the problem is warm-started. Scaling the problem is useful to prevent 
 * large changes in the active set and to guard against ill-conditioning in the objective function. 
 *  
 */


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include <stdio.h>
#include <qpalm/scaling.h>
#include <qpalm/lin_alg.h>
#include <ladel.h>

// Set values lower than threshold MIN_SCALING to 1
void limit_scaling(c_float *D, size_t n) 
{
    size_t i;

    for (i = 0; i < n; i++) {
        D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
    }
}

void scale_data(QPALMWorkspace *work) 
{
    size_t n = work->data->n;
    size_t m = work->data->m;
    vec_set_scalar(work->scaling->D, 1, n);
    vec_set_scalar(work->scaling->E, 1, m);
    
    c_int i;
    //Ruiz on constraint matrix A
    for (i = 0; i < work->settings->scaling; i++) {

        // Set D_temp = vecnorm(A,inf,1) (cols) and E_temp = vecnorm(A,inf,2) (rows)
        mat_inf_norm_cols(work->data->A, work->D_temp);
        mat_inf_norm_rows(work->data->A, work->E_temp);

        // Set to 1 values with 0 norms (avoid crazy scaling)
        limit_scaling(work->D_temp, n);
        limit_scaling(work->E_temp, m);

        // Take square root of norms
        vec_ew_sqrt(work->D_temp, work->D_temp, n);
        vec_ew_sqrt(work->E_temp, work->E_temp, m);

        // 1./D and 1./E
        vec_ew_recipr(work->D_temp, work->D_temp, n);
        vec_ew_recipr(work->E_temp, work->E_temp, m);

        // Equilibrate matrix A
        // A <- EAD
        
        ladel_scale_rows(work->data->A, work->solver->E_temp);
        ladel_scale_columns(work->data->A, work->solver->D_temp);
        // Update equilibration matrices D and E
        vec_ew_prod(work->scaling->D, work->D_temp, work->scaling->D, n);
        vec_ew_prod(work->scaling->E, work->E_temp, work->scaling->E, m);
    }

    // Equilibrate matrix Q and vector q
    // Q <- cDQD, q <- cDq
    vec_ew_prod(work->scaling->D, work->data->q, work->data->q, n);
    vec_ew_prod(work->scaling->D, work->Qx, work->Qx, n);
    prea_vec_copy(work->scaling->D, work->D_temp, n);
    vec_add_scaled(work->Qx, work->data->q, work->dphi, 1, n);
    work->scaling->c = 1/c_max(1.0, vec_norm_inf(work->dphi, n));
    vec_self_mult_scalar(work->data->q, work->scaling->c, n);
    vec_self_mult_scalar(work->Qx, work->scaling->c, n);
    
    ladel_scale_columns(work->data->Q, work->solver->D_temp);
    ladel_scale_rows(work->data->Q, work->solver->D_temp);
    ladel_scale_scalar(work->data->Q, work->scaling->c);

    // Store cinv, Dinv, Einv
    vec_ew_recipr(work->scaling->D, work->scaling->Dinv, n);
    vec_ew_recipr(work->scaling->E, work->scaling->Einv, m);
    work->scaling->cinv = (c_float) 1.0/work->scaling->c;

    // Scale problem vectors l, u
    vec_ew_prod(work->scaling->E, work->data->bmin, work->data->bmin, m);
    vec_ew_prod(work->scaling->E, work->data->bmax, work->data->bmax, m);    

    // Scale initial vectors x, Ax and y (Qx already scaled)
    vec_ew_prod(work->x, work->scaling->Dinv, work->x, n);
    vec_ew_prod(work->Ax, work->scaling->E, work->Ax, m);
    vec_ew_prod(work->y, work->scaling->Einv, work->y, m);
    vec_self_mult_scalar(work->y, work->scaling->c, m);
}

void unscale_data(QPALMWorkspace *work)
{
    size_t n = work->data->n;
    size_t m = work->data->m;
    if (work->settings->scaling) 
    {
        // A (was scaled to EAD)
        ladel_scale_rows(work->data->A, work->scaling->Einv);
        ladel_scale_columns(work->data->A, work->scaling->Dinv);

        // Q (was scaled to cDQD)
        ladel_scale_columns(work->data->Q, work->scaling->Dinv);
        ladel_scale_rows(work->data->Q, work->scaling->Dinv);
        ladel_scale_scalar(work->data->Q, work->scaling->cinv);

        // q (was scaled to cDq)
        vec_ew_prod(work->scaling->Dinv, work->data->q, work->data->q, n);
        vec_self_mult_scalar(work->data->q, work->scaling->cinv, n);

        // bmin/bmax (was scaled to Ebmin/Ebmax)
        vec_ew_prod(work->scaling->Einv, work->data->bmin, work->data->bmin, m);
        vec_ew_prod(work->scaling->Einv, work->data->bmax, work->data->bmax, m);    
    }
}


# ifdef __cplusplus
}
# endif // ifdef __cplusplus
