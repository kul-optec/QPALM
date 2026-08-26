/**
 * @file iteration.c
 * @author Ben Hermans
 * @brief QPALM main solver routines.
 * @details This file contains the functions that make up the qpalm algorithm (the functions that qpalm_solve will use). 
 * These include the computation of the residuals at the start of the iteration, 
 * the update of the primal variables in an inner iteration, 
 * the update of the penalty factors and dual variables in an outer iteration,
 * the computation of the primal and dual objective values, etc.
 */

#include <qpalm/iteration.h>
#include <qpalm/lin_alg.h>
#include <qpalm/solver_interface.h>
#include <qpalm/newton.h>
#include <qpalm/linesearch.h>
#include <qpalm/nonconvex.h>
#include <qpalm/util.h>

#include <ladel.h>

void compute_residuals(QPALMWorkspace *work, solver_common *c) {

    //Axys = Ax + y./sigma
    vec_ew_prod(work->y, work->sigma_inv, work->temp_m, work->data->m);
    vec_add_scaled(work->Ax, work->temp_m, work->Axys, 1, work->data->m);
    //z = min(max(Axys,bmin),bmax)
    vec_ew_mid_vec(work->Axys, work->data->bmin, work->data->bmax, work->z, work->data->m);
    //pri_res = Ax-z
    vec_add_scaled(work->Ax, work->z, work->pri_res, -1, work->data->m);
    //yh = y + pri_res.*sigma
    vec_ew_prod(work->pri_res, work->sigma, work->temp_m, work->data->m);
    vec_add_scaled(work->y, work->temp_m, work->yh, 1, work->data->m);
    //df = Qx + q
    vec_add_scaled(work->Qx, work->data->q, work->df, 1, work->data->n);
    
    if (work->settings->proximal) {
      //df = Qx + q +1/gamma*(x-x0)
      // NB work->Qx contains Qx+1/gamma*x
      vec_add_scaled(work->df, work->x0, work->df, -1/work->gamma, work->data->n);
    }
    // Atyh = A'*yh
    mat_tpose_vec(work->data->A, work->solver->yh, work->solver->Atyh, c);
    //dphi = df+Atyh
    vec_add_scaled(work->df, work->Atyh, work->dphi, 1, work->data->n);
}

void initialize_sigma(QPALMWorkspace *work, solver_common *c) {

    // Compute initial sigma
    size_t n = work->data->n;
    size_t m = work->data->m;
    c_float f = 0.5*vec_prod(work->x, work->Qx, n) + vec_prod(work->data->q, work->x, n);
    if (work->settings->proximal) {
        f -= 0.5/work->gamma*vec_prod(work->x, work->x, n);
    }
    vec_ew_mid_vec(work->Ax, work->data->bmin, work->data->bmax, work->temp_m, m);
    vec_add_scaled(work->Ax, work->temp_m, work->temp_m, -1, m);
    c_float dist2 = vec_prod(work->temp_m, work->temp_m, m);
    vec_set_scalar(work->sigma, c_max(1e-4, c_min(work->settings->sigma_init*c_max(1,c_absval(f))/c_max(1,0.5*dist2),1e4)), m);

    // Set fields related to sigma
    vec_ew_recipr(work->sigma, work->sigma_inv, m);
    vec_ew_sqrt(work->sigma, work->sqrt_sigma, m);
    work->sqrt_sigma_max = c_sqrt(work->settings->sigma_max);
    
    if (work->solver->factorization_method == FACTORIZE_SCHUR)
    {
        work->solver->At_sqrt_sigma = ladel_sparse_free(work->solver->At_sqrt_sigma);
        work->solver->At_sqrt_sigma = ladel_transpose(work->data->A, TRUE, c);
        ladel_scale_columns(work->solver->At_sqrt_sigma, work->sqrt_sigma);
    }

}

void update_sigma(QPALMWorkspace* work, solver_common *c) {
    
    work->nb_sigma_changed = 0;
    c_float *At_scalex = work->solver->At_scale;
    c_float pri_res_unscaled_norm = vec_norm_inf(work->pri_res, work->data->m);
    c_float sigma_temp, mult_factor;
    c_int *sigma_changed = work->solver->enter;
    size_t k;
    for (k = 0; k < work->data->m; k++) {
        if ((c_absval(work->pri_res[k]) > work->settings->theta*c_absval(work->pri_res_in[k])) && work->solver->active_constraints[k]) {
            mult_factor = c_max(1.0, work->settings->delta * c_absval(work->pri_res[k]) / (pri_res_unscaled_norm + 1e-6));
            sigma_temp = mult_factor * work->sigma[k];
            if (sigma_temp <= work->settings->sigma_max) { 
                if (work->sigma[k] != sigma_temp) {
                    sigma_changed[work->nb_sigma_changed] = (c_int)k;
                    work->nb_sigma_changed++;
                }               
                work->sigma[k] = sigma_temp;
                work->sigma_inv[k] = 1.0/sigma_temp;
                mult_factor = c_sqrt(mult_factor);
                work->sqrt_sigma[k] = mult_factor * work->sqrt_sigma[k];
                At_scalex[k] = mult_factor;
            } else {
                if (work->sigma[k] != work->settings->sigma_max) {
                    sigma_changed[work->nb_sigma_changed] = (c_int)k;
                    work->nb_sigma_changed++;
                } 
                work->sigma[k] = work->settings->sigma_max;
                work->sigma_inv[k] = 1.0/work->settings->sigma_max;
                At_scalex[k] = work->sqrt_sigma_max / work->sqrt_sigma[k];
                work->sqrt_sigma[k] = work->sqrt_sigma_max;
            }
        } else {
            At_scalex[k] = 1.0;
        }
    }

    // TODO implement updating sigma in KKT system
    if (work->solver->factorization_method == FACTORIZE_SCHUR)
        ladel_scale_columns(work->solver->At_sqrt_sigma, work->solver->At_scale);

    if (work->solver->first_factorization || (work->settings->proximal && work->gamma < work->settings->gamma_max) || 
        (work->nb_sigma_changed > 
            c_min(work->settings->max_rank_update_fraction*(work->data->n+work->data->m), 0.25*work->settings->max_rank_update))) 
    {
        work->solver->reset_newton = TRUE;
    } else if (work->nb_sigma_changed == 0){
        /* do nothing */
    } else {
        ldlupdate_sigma_changed(work, c);
    }
}

void update_gamma(QPALMWorkspace *work) {
    
    if (work->gamma < work->settings->gamma_max) 
    {
        c_float prev_gamma = work->gamma;
        work->gamma = c_min(work->gamma*work->settings->gamma_upd, work->settings->gamma_max);
        work->solver->reset_newton = TRUE;
        vec_add_scaled(work->Qx, work->x, work->Qx, 1/work->gamma - 1/prev_gamma, work->data->n);
    }
    
}

void boost_gamma(QPALMWorkspace *work, solver_common *c) {

    c_float prev_gamma = work->gamma;
    if (work->solver->nb_active_constraints) {
        // work->gamma = 1e10;
        solver_sparse *AtsigmaA = NULL;
        size_t nb_active = 0;
        for (size_t i = 0; i < work->data->m; i++){
            if (work->solver->active_constraints[i]){
                work->solver->enter[nb_active] = (c_int)i;
                nb_active++;
            }      
        }
        solver_sparse *A = NULL, *At = NULL;
        if (work->solver->factorization_method == FACTORIZE_KKT)
        {
            work->gamma = 1e10;
            // size_t i;
            // At = ladel_column_submatrix(work->solver->At, work->solver->enter, nb_active);
            // A = ladel_transpose(At, TRUE, c);
            // for (i = 0; i < nb_active; i++) 
            //     work->temp_m[i] = work->sigma[work->solver->enter[i]];
            // AtsigmaA = ladel_mat_diag_mat_transpose(At, A, work->temp_m, c);
        } else if (work->solver->factorization_method == FACTORIZE_SCHUR)
        {
            At = ladel_column_submatrix(work->solver->At_sqrt_sigma, work->solver->enter, nb_active);
            A = ladel_transpose(At, TRUE, c);
            AtsigmaA = ladel_mat_mat_transpose(At, A, c);
            work->gamma = c_max(work->settings->gamma_max, 1e14/gershgorin_max(AtsigmaA, work->temp_n, work->neg_dphi));
        }

        // work->gamma = c_max(work->settings->gamma_max, 1e14/gershgorin_max(AtsigmaA, work->temp_n, work->neg_dphi));
        work->gamma_maxed = TRUE;
        A = ladel_sparse_free(A);
        At = ladel_sparse_free(At);
        AtsigmaA = ladel_sparse_free(AtsigmaA);
    } else {
        work->gamma = 1e12;
    }
    if (prev_gamma != work->gamma) {
        vec_add_scaled(work->Qx, work->x, work->Qx, 1.0/work->gamma - 1.0/prev_gamma, work->data->n);
        vec_add_scaled(work->Qd, work->d, work->Qd, work->tau/work->gamma - work->tau/prev_gamma, work->data->n);
        work->solver->reset_newton = TRUE;
    }
}

void update_or_boost_gamma(QPALMWorkspace *work, solver_common *c, c_int iter_out)
{
    if (!work->gamma_maxed && iter_out > 0 && work->solver->nb_enter == 0 
        && work->solver->nb_leave == 0 && work->info->pri_res_norm < work->eps_pri) 
    {
        //Axys = Ax + y./sigma
        vec_ew_div(work->y, work->sigma, work->temp_m, work->data->m);
        vec_add_scaled(work->Ax, work->temp_m, work->Axys, 1, work->data->m);
        set_active_constraints(work);
        set_entering_leaving_constraints(work);
        if (work->solver->nb_enter == 0 && work->solver->nb_leave == 0) 
        {
            boost_gamma(work, c);
        } 
        else 
        {
            update_gamma(work);
        }
    } 
    else 
    {
        update_gamma(work);
    }
}

void update_proximal_point_and_penalty(QPALMWorkspace *work, solver_common *c, c_int iter_out, c_float *eps_k_abs, c_float *eps_k_rel)
{
    if (work->settings->nonconvex)
    {
        c_float eps_k; 
        c_int m = work->data->m;
        if (work->settings->scaling) 
        {
            /**NB Implementation detail: store Einv*Ax and Einv*z in temp_2m. 
             * The infinity norm of that vector is then equal to the maximum of both norms. */
            vec_ew_prod(work->scaling->Einv, work->Ax, work->temp_2m, m);
            vec_ew_prod(work->scaling->Einv, work->z, work->temp_2m + m, m);
            eps_k =  (*eps_k_abs) + (*eps_k_rel)*vec_norm_inf(work->temp_2m, 2 * m);
        } 
        else 
        {
            eps_k =  (*eps_k_abs) + (*eps_k_rel)*c_max(vec_norm_inf(work->Ax, m), vec_norm_inf(work->z, m));
        }

        if (work->info->pri_res_norm < eps_k)
        {
            prea_vec_copy(work->x, work->x0, work->data->n);
            *eps_k_abs = c_max(work->settings->eps_abs, work->settings->rho*(*eps_k_abs));
            *eps_k_rel = c_max(work->settings->eps_rel, work->settings->rho*(*eps_k_rel));
        }
        else
        {
            /* We do not update the tolerances and proximal point in this case */
        }
    } 
    else
    {
        if(work->settings->proximal) 
        {
            update_or_boost_gamma(work, c, iter_out);                
            prea_vec_copy(work->x, work->x0, work->data->n);
        }
    }
}

void update_dual_iterate_and_parameters(QPALMWorkspace *work, solver_common *c, c_int iter_out, c_float *eps_k_abs, c_float *eps_k_rel)
{
    c_int n = work->data->n, m = work->data->m;

    if (iter_out > 0 && work->info->pri_res_norm > work->eps_pri) 
    {
        update_sigma(work, c);
    }

    prea_vec_copy(work->yh, work->y, m);
    prea_vec_copy(work->Atyh, work->Aty, n);

    work->eps_abs_in = c_max(work->settings->eps_abs, work->settings->rho*work->eps_abs_in);
    work->eps_rel_in = c_max(work->settings->eps_rel, work->settings->rho*work->eps_rel_in); 

    update_proximal_point_and_penalty(work, c, iter_out, eps_k_abs, eps_k_rel);

    prea_vec_copy(work->pri_res, work->pri_res_in, m);
}

void update_primal_iterate(QPALMWorkspace *work, solver_common *c) {
    
    newton_set_direction(work, c);

    work->tau = exact_linesearch(work, c);

    //x_prev = x
    prea_vec_copy(work->x, work->x_prev, work->data->n);
    //dphi_prev = dphi 
    prea_vec_copy(work->dphi, work->dphi_prev, work->data->n);
    //x = x+tau*d
    vec_add_scaled(work->x, work->d, work->x, work->tau, work->data->n);
    vec_self_mult_scalar(work->Qd, work->tau, work->data->n); //Qdx used in dua_infeas check
    vec_self_mult_scalar(work->Ad, work->tau, work->data->m); //Adx used in dua_infeas check
    vec_add_scaled(work->Qx, work->Qd, work->Qx, 1, work->data->n);
    vec_add_scaled(work->Ax, work->Ad, work->Ax, 1, work->data->m);
}

c_float compute_objective(QPALMWorkspace *work) {
    
    c_float objective = 0;
    size_t n = work->data->n;
    size_t i = 0; 

    if (work->settings->proximal) {
        if(n >= 4) {
            for (; i <= n-4; i+=4) {
                objective +=  (0.5*(work->Qx[i] - 1/work->gamma*work->x[i]) + work->data->q[i])*work->x[i] 
                            + (0.5*(work->Qx[i+1] - 1/work->gamma*work->x[i+1]) + work->data->q[i+1])*work->x[i+1]
                            + (0.5*(work->Qx[i+2] - 1/work->gamma*work->x[i+2]) + work->data->q[i+2])*work->x[i+2] 
                            + (0.5*(work->Qx[i+3] - 1/work->gamma*work->x[i+3]) + work->data->q[i+3])*work->x[i+3];
            }
        }
        for (; i < n; i++) {
            objective += (0.5*(work->Qx[i] - 1/work->gamma*work->x[i])+ work->data->q[i])*work->x[i];
        }
    } else {
        if(n >= 4) {
            for (; i <= n-4; i+=4) {
                objective += (0.5*work->Qx[i] + work->data->q[i])*work->x[i] 
                            + (0.5*work->Qx[i+1] + work->data->q[i+1])*work->x[i+1]
                            + (0.5*work->Qx[i+2] + work->data->q[i+2])*work->x[i+2] 
                            + (0.5*work->Qx[i+3] + work->data->q[i+3])*work->x[i+3];
            }
        }
        for (; i < n; i++) {
            objective += (0.5*work->Qx[i] + work->data->q[i])*work->x[i];
        }
    }

    if (work->settings->scaling) {
        objective *= work->scaling->cinv;
    }

    objective += work->data->c;

    return objective;
}

c_float compute_dual_objective(QPALMWorkspace *work, solver_common *c) {

    c_float dual_objective = 0;

    vec_add_scaled(work->Atyh, work->data->q, work->neg_dphi, 1.0, work->data->n);
    ladel_dense_solve(work->solver->LD_Q, work->neg_dphi, work->D_temp, c);

    dual_objective -= 0.5*vec_prod(work->neg_dphi, work->D_temp, work->data->n);
    for (size_t i = 0; i < work->data->m; i++) {
      dual_objective -= work->yh[i] > 0 ? work->yh[i]*work->data->bmax[i] : work->yh[i]*work->data->bmin[i];
    }

    if(work->settings->scaling) {
      dual_objective *= work->scaling->cinv;
    }

    dual_objective += work->data->c;

    return dual_objective;
}
