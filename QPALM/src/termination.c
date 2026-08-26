/**
 * @file termination.c
 * @author Ben Hermans
 * @brief Routines to check the termination and infeasibility criteria.
 * @details The routines in this file compute the primal and dual residuals, 
 * the primal and dual tolerances, check whether the problem is solved 
 * completely, unscale and store the solution if that is the case, check
 *  whether the intermediate problem is solved and whether one of the 
 * infeasibility criteria hold. In other words, all routines related to 
 * the termination of the optimization algorithm are grouped in this file.
 */
#include <qpalm/termination.h>
#include <qpalm/lin_alg.h>
#include <qpalm/constants.h>
#include <qpalm/global_opts.h>
#include <qpalm/util.h>
#include <qpalm/iteration.h>

void calculate_residual_norms_and_tolerances(QPALMWorkspace *work) {
    calculate_primal_residual(work);
    calculate_dual_residuals(work);
    calculate_primal_tolerance(work);
    calculate_dual_tolerances(work);
}

void calculate_primal_residual(QPALMWorkspace *work) {
    size_t m = work->data->m;
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->pri_res, work->temp_m, m);
        work->info->pri_res_norm = vec_norm_inf(work->temp_m, m);
    } else {
        work->info->pri_res_norm = vec_norm_inf(work->pri_res, m);
    }
}

void calculate_dual_residuals(QPALMWorkspace *work) {
    size_t n = work->data->n;
    if (work->settings->scaling) {
        if (work->settings->proximal) {
            vec_add_scaled(work->x, work->x0, work->xx0, -1, n);
            vec_add_scaled(work->dphi, work->xx0, work->temp_n, -1/work->gamma, n);
            vec_ew_prod(work->scaling->Dinv, work->temp_n, work->temp_n, n);
            work->info->dua_res_norm = vec_norm_inf(work->temp_n, n);
            vec_ew_prod(work->scaling->Dinv, work->dphi, work->temp_n, n);
            work->info->dua2_res_norm = vec_norm_inf(work->temp_n, n);
        } else {
            vec_ew_prod(work->scaling->Dinv, work->dphi, work->temp_n, n);
            work->info->dua_res_norm = vec_norm_inf(work->temp_n, n);
            work->info->dua2_res_norm = work->info->dua_res_norm;
        } 
        work->info->dua_res_norm *= work->scaling->cinv;
        work->info->dua2_res_norm *= work->scaling->cinv;
            
    } else {
        if (work->settings->proximal) {
            vec_add_scaled(work->x, work->x0, work->xx0, -1, n);
            vec_add_scaled(work->dphi, work->xx0, work->temp_n, -1/work->gamma, n);
            work->info->dua_res_norm = vec_norm_inf(work->temp_n, n);
            work->info->dua2_res_norm = vec_norm_inf(work->dphi, n);
        } else {
            work->info->dua_res_norm = vec_norm_inf(work->dphi, n);
            work->info->dua2_res_norm = work->info->dua_res_norm;
        }
    }
}

void calculate_primal_tolerance(QPALMWorkspace *work) {
    size_t m = work->data->m;
    if (work->settings->scaling) {
        /**NB Implementation detail: store Einv*Ax and Einv*z in temp_2m. 
         * The infinity norm of that vector is equal to the maximum
         * of the infinity norms of Einv*Ax and Einv*z.*/
        vec_ew_prod(work->scaling->Einv, work->Ax, work->temp_2m, m);
        vec_ew_prod(work->scaling->Einv, work->z, work->temp_2m + m, m);
        work->eps_pri =  work->settings->eps_abs + work->settings->eps_rel*vec_norm_inf(work->temp_2m, 2 * m);
    } else {
        work->eps_pri =  work->settings->eps_abs + work->settings->eps_rel*c_max(
                                vec_norm_inf(work->Ax, m),
                                vec_norm_inf(work->z, m));
    }
}

void calculate_dual_tolerances(QPALMWorkspace *work) {
    size_t n = work->data->n;
    c_float norm_DinvQx, norm_Dinvq, norm_DinvAtyh, max_norm;
    if (work->settings->scaling) {
        if (work->settings->proximal) {
            vec_add_scaled(work->Qx, work->x, work->temp_n, -1/work->gamma, n);
            vec_ew_prod(work->scaling->Dinv, work->temp_n, work->temp_n, n);
        } else {
            vec_ew_prod(work->scaling->Dinv, work->Qx, work->temp_n, n);
        }
        norm_DinvQx = vec_norm_inf(work->temp_n, n);
        vec_ew_prod(work->scaling->Dinv, work->data->q, work->temp_n, n);
        norm_Dinvq = vec_norm_inf(work->temp_n, n);
        vec_ew_prod(work->scaling->Dinv, work->Atyh, work->temp_n, n);
        norm_DinvAtyh = vec_norm_inf(work->temp_n, n);
    } else {
        if (work->settings->proximal) {
            vec_add_scaled(work->Qx, work->x, work->temp_n, -1/work->gamma, n);
            norm_DinvQx = vec_norm_inf(work->temp_n, n);
        } else {
            norm_DinvQx = vec_norm_inf(work->Qx, n);
        }
        norm_Dinvq = vec_norm_inf(work->data->q, n);
        norm_DinvAtyh = vec_norm_inf(work->Atyh, n);
    }
    
    max_norm = c_max(norm_DinvQx, c_max(norm_Dinvq, norm_DinvAtyh));
    if (work->settings->scaling) max_norm *= work->scaling->cinv;

    work->eps_dua = work->settings->eps_abs + work->settings->eps_rel*max_norm;
    work->eps_dua_in = work->eps_abs_in + work->eps_rel_in*max_norm;
}

c_int is_solved(QPALMWorkspace *work) {
    return (work->info->pri_res_norm < work->eps_pri) && 
            (work->info->dua_res_norm < work->eps_dua);
}

c_int is_primal_infeasible(QPALMWorkspace *work) {
    size_t n = work->data->n;
    size_t m = work->data->m;
    c_float eps_pinf_norm_Edy;

    //dy = yh-y
    vec_add_scaled(work->yh, work->y, work->delta_y, -1, m);
    if (work->settings->scaling) {
        //Edy = E.*dy
        vec_ew_prod(work->scaling->E, work->delta_y, work->temp_m, m);
        eps_pinf_norm_Edy = work->settings->eps_prim_inf*vec_norm_inf(work->temp_m, m);
    } else {
        eps_pinf_norm_Edy = work->settings->eps_prim_inf*vec_norm_inf(work->delta_y, m);
    }
     
    if (eps_pinf_norm_Edy == 0) { //dy == 0
        return 0;
    }

    vec_add_scaled(work->Atyh, work->Aty, work->Atdelta_y, -1, n);
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->Atdelta_y, work->Atdelta_y, n);
    }

    //out_of_bounds = bmax'*max(dy,0) + bmin'*min(dy,0)
    c_float out_of_bounds = 0;
    if (work->settings->scaling)
    {
        for(size_t i=0; i < m; i++) 
        {
            out_of_bounds += (work->data->bmax[i] < work->scaling->E[i]*QPALM_INFTY) ? work->data->bmax[i]*c_max(work->delta_y[i], 0) : 0;
            out_of_bounds += (work->data->bmin[i] > -work->scaling->E[i]*QPALM_INFTY) ? work->data->bmin[i]*c_min(work->delta_y[i], 0) : 0;
        }
    }
    else
    {
        for(size_t i=0; i < m; i++) 
        {
            out_of_bounds += (work->data->bmax[i] < QPALM_INFTY) ? work->data->bmax[i]*c_max(work->delta_y[i], 0) : 0;
            out_of_bounds += (work->data->bmin[i] > -QPALM_INFTY) ? work->data->bmin[i]*c_min(work->delta_y[i], 0) : 0;
        }
    }

    return (vec_norm_inf(work->Atdelta_y, n) <= eps_pinf_norm_Edy)
        && (out_of_bounds <= -eps_pinf_norm_Edy);
    
}

c_int is_dual_infeasible(QPALMWorkspace *work) {
    c_float eps_dinf_norm_Ddx, dxQdx, dxdx;
    size_t n = work->data->n;
    size_t m = work->data->m;

    //dx = x-x_prev
    vec_add_scaled(work->x, work->x_prev, work->delta_x, -1, n);
    if (work->settings->scaling) {
        //D*dx
        vec_ew_prod(work->scaling->D, work->delta_x, work->temp_n, n);
        eps_dinf_norm_Ddx = work->settings->eps_dual_inf*vec_norm_inf(work->temp_n, n);
        dxdx = vec_prod(work->temp_n, work->temp_n, n);
    } else {
        eps_dinf_norm_Ddx = work->settings->eps_dual_inf*vec_norm_inf(work->delta_x, n);
        dxdx = vec_prod(work->delta_x, work->delta_x, n);
    }
    
    if (eps_dinf_norm_Ddx == 0) { //dx == 0
        return 0;
    }

    size_t k;
    //NB Adx = work->Ad (= tau*Ad of the previous iteration)
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->Ad, work->Adelta_x, m);
        for (k = 0; k < m; k++) {
            if ((work->data->bmax[k] < work->scaling->E[k]*QPALM_INFTY && work->Adelta_x[k] >= eps_dinf_norm_Ddx)
                || (work->data->bmin[k] > -work->scaling->E[k]*QPALM_INFTY && work->Adelta_x[k] <= -eps_dinf_norm_Ddx)) {
                return 0;
            }
        }      
    } else {
        for (k = 0; k < m; k++) {
            if ((work->data->bmax[k] < QPALM_INFTY && work->Ad[k] >= eps_dinf_norm_Ddx)
                || (work->data->bmin[k] > -QPALM_INFTY && work->Ad[k] <= -eps_dinf_norm_Ddx)) {
                return 0;
            }
        }
    }
    //NB Qdx = work->Qd (= tau*Qd of the previous iteration)
    //NB Qdx = work->Qd - tau/gamma*d (= tau*Qd of the previous iteration) if proximal is used
    if (work->settings->proximal) {
        vec_add_scaled(work->Qd, work->d, work->temp_n, -work->tau/work->gamma, n);
        dxQdx = vec_prod(work->delta_x, work->temp_n, n);
    } else {
        dxQdx = vec_prod(work->Qd, work->delta_x, n);
    }
    if (work->settings->scaling) {
        return (dxQdx <= -work->scaling->c*work->settings->eps_dual_inf*work->settings->eps_dual_inf*dxdx)
            || ((dxQdx <= work->scaling->c*work->settings->eps_dual_inf*work->settings->eps_dual_inf*dxdx)
                && (vec_prod(work->data->q, work->delta_x, n) <= -work->scaling->c*eps_dinf_norm_Ddx));
    } else {
        return (dxQdx <= -work->settings->eps_dual_inf*work->settings->eps_dual_inf*dxdx)
            || ((dxQdx <= work->settings->eps_dual_inf*work->settings->eps_dual_inf*dxdx)
                && (vec_prod(work->data->q, work->delta_x, n) <= -eps_dinf_norm_Ddx));
    }  
}

void store_solution(QPALMWorkspace *work) {
    if (work->settings->scaling) {
        vec_ew_prod(work->x, work->scaling->D, work->solution->x, work->data->n);
        vec_ew_prod(work->yh, work->scaling->E, work->solution->y, work->data->m);
        vec_self_mult_scalar(work->solution->y, work->scaling->cinv, work->data->m);
    } else {
        prea_vec_copy(work->x, work->solution->x, work->data->n);
        prea_vec_copy(work->yh, work->solution->y, work->data->m);
    }
    work->info->objective = compute_objective(work);
}

c_int check_subproblem_termination(QPALMWorkspace *work) {
    return (work->info->dua2_res_norm <= work->eps_dua_in);
}
