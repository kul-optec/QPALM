/**
 * @file qpalm.c
 * @author Ben Hermans
 * @brief QPALM main solver API.
 * @details This file contains the main functions that can be called by the user.
 * The user can load the default settings, setup the workspace with data and settings,
 * run the solver, and cleanup the workspace afterwards.
 */
# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include <qpalm.h>
#include <qpalm/global_opts.h>
#include <qpalm/constants.h>
#include <qpalm/validate.h>
#include <qpalm/lin_alg.h>
#include <qpalm/util.h>
#include <qpalm/scaling.h>
#include <qpalm/linesearch.h>
#include <qpalm/termination.h>
#include <qpalm/solver_interface.h>
#include <qpalm/newton.h>
#include <qpalm/nonconvex.h>
#include <qpalm/iteration.h>

#include <ladel.h>

/**********************
* Main API Functions *
**********************/

void qpalm_set_default_settings(QPALMSettings *settings) 
{
    settings->max_iter                  = MAX_ITER;                           /* maximum iterations */
    settings->inner_max_iter            = INNER_MAX_ITER;                     /* maximum iterations per subproblem */
    settings->eps_abs                   = (c_float)EPS_ABS;                   /* absolute convergence tolerance */
    settings->eps_rel                   = (c_float)EPS_REL;                   /* relative convergence tolerance */
    settings->eps_abs_in                = (c_float)EPS_ABS_IN;                /* intermediate absolute convergence tolerance */
    settings->eps_rel_in                = (c_float)EPS_REL_IN;                /* intermediate relative convergence tolerance */
    settings->rho                       = (c_float)RHO;                       /* tolerance scaling factor */
    settings->eps_prim_inf              = (c_float)EPS_PRIM_INF;              /* primal infeasibility tolerance */
    settings->eps_dual_inf              = (c_float)EPS_DUAL_INF;              /* dual infeasibility tolerance */
    settings->theta                     = (c_float)THETA;                     /* penalty update criterion parameter */
    settings->delta                     = (c_float)DELTA;                     /* penalty update factor */
    settings->sigma_max                 = (c_float)SIGMA_MAX;                 /* penalty parameter cap */
    settings->sigma_init                = (c_float)SIGMA_INIT;                /* initial penalty parameter (guideline) */
    settings->proximal                  = PROXIMAL;                           /* boolean, proximal method of multipliers*/
    settings->gamma_init                = (c_float)GAMMA_INIT;                /* proximal penalty parameter */
    settings->gamma_upd                 = (c_float)GAMMA_UPD;                 /* proximal penalty update factor*/
    settings->gamma_max                 = (c_float)GAMMA_MAX;                 /* proximal penalty parameter cap*/
    settings->scaling                   = SCALING;                            /* boolean, scaling */
    settings->nonconvex                 = NONCONVEX;                          /* boolean, nonconvex */
    settings->warm_start                = WARM_START;                         /* boolean, warm start solver */
    settings->verbose                   = VERBOSE;                            /* boolean, write out progress */
    settings->print_iter                = PRINT_ITER;                         /* frequency of printing */
    settings->reset_newton_iter         = RESET_NEWTON_ITER;                  /* frequency of performing a full Cholesky factorization */
    settings->enable_dual_termination   = ENABLE_DUAL_TERMINATION;            /* allow for dual termination (useful in branch and bound) */
    settings->dual_objective_limit      = (c_float)DUAL_OBJECTIVE_LIMIT;      /* termination value for the dual objective (useful in branch and bound) */
    settings->time_limit                = (c_float)TIME_LIMIT;                /* time limit */
    settings->ordering                  = ORDERING;                           /* ordering */
    settings->factorization_method      = FACTORIZATION_METHOD;               /* factorization method (kkt or schur) */
    settings->max_rank_update           = MAX_RANK_UPDATE;                    /* maximum rank of the update */
    settings->max_rank_update_fraction  = (c_float)MAX_RANK_UPDATE_FRACTION;  /* maximum rank (relative to n+m) of the update */
}


QPALMWorkspace* qpalm_setup(const QPALMData *data, const QPALMSettings *settings) 
{
    QPALMWorkspace *work; // Workspace

    // Validate data
    if (!validate_data(data)) 
    {
        # ifdef QPALM_PRINTING
        qpalm_eprint("Data validation returned failure");
        # endif /* ifdef QPALM_PRINTING */
        return QPALM_NULL;
    }

    // Validate settings
    if (!validate_settings(settings)) 
    {
        # ifdef QPALM_PRINTING
        qpalm_eprint("Settings validation returned failure");
        # endif /* ifdef QPALM_PRINTING */
        return QPALM_NULL;
    }

    // Allocate empty workspace
    work = qpalm_calloc(1, sizeof(QPALMWorkspace));

    if (!work) 
    {
        # ifdef QPALM_PRINTING
        qpalm_eprint("allocating work failure");
        # endif /* ifdef QPALM_PRINTING */
        return QPALM_NULL;
    }

    // Start and allocate directly timer
    # ifdef QPALM_TIMING
    work->timer = qpalm_malloc(sizeof(QPALMTimer));
    qpalm_tic(work->timer);
    # endif /* ifdef QPALM_TIMING */

    // Copy settings
    work->settings = copy_settings(settings);
    work->sqrt_delta = c_sqrt(work->settings->delta);
    work->gamma = work->settings->gamma_init;

    size_t n = data->n;
    size_t m = data->m;

    //Initialize the solver for the linear system
    work->solver = qpalm_calloc(1, sizeof(QPALMSolver));
    solver_common common, *c;
    c = &common;

    // Copy problem data into workspace
    work->data       = qpalm_calloc(1, sizeof(QPALMData));
    work->data->n    = data->n;           
    work->data->m    = data->m;                   
    work->data->bmin = vec_copy(data->bmin, m);      
    work->data->bmax = vec_copy(data->bmax, m);       
    work->data->q    = vec_copy(data->q, n);          
    work->data->c    = data->c;

    work->data->A = ladel_sparse_allocate_and_copy(data->A); 
    work->data->Q = ladel_sparse_allocate_and_copy(data->Q);
    ladel_to_upper_diag(work->data->Q);

    // Allocate internal solver variables 
    work->x        = qpalm_calloc(n, sizeof(c_float));
    work->y        = qpalm_calloc(m, sizeof(c_float));
    work->Ax       = qpalm_calloc(m, sizeof(c_float));
    work->Qx       = qpalm_calloc(n, sizeof(c_float));
    work->x_prev   = qpalm_calloc(n, sizeof(c_float));
    work->Aty      = qpalm_calloc(n, sizeof(c_float));

    work->x0 = qpalm_calloc(n, sizeof(c_float));

    work->initialized = FALSE;

    // Workspace variables
    work->temp_m   = qpalm_calloc(m, sizeof(c_float));
    work->temp_n   = qpalm_calloc(n, sizeof(c_float));
    work->sigma = qpalm_calloc(m, sizeof(c_float));
    work->sigma_inv = qpalm_calloc(m, sizeof(c_float));
    work->nb_sigma_changed = 0;

    work->z  = qpalm_calloc(m, sizeof(c_float));
    work->Axys = qpalm_calloc(m, sizeof(c_float));
    work->pri_res = qpalm_calloc(m, sizeof(c_float));
    work->pri_res_in = qpalm_calloc(m, sizeof(c_float));
    work->df = qpalm_calloc(n, sizeof(c_float));
    
    work->xx0 = qpalm_calloc(n, sizeof(c_float));
    work->dphi = qpalm_calloc(n, sizeof(c_float));
    work->dphi_prev = qpalm_calloc(n, sizeof(c_float));

    // Linesearch variables
    work->sqrt_sigma  = qpalm_calloc(m, sizeof(c_float));
    work->delta       = qpalm_calloc(m*2, sizeof(c_float));
    work->alpha       = qpalm_calloc(m*2, sizeof(c_float));
    work->delta2      = qpalm_calloc(m*2, sizeof(c_float));
    work->delta_alpha = qpalm_calloc(m*2, sizeof(c_float));
    work->temp_2m     = qpalm_calloc(m*2, sizeof(c_float));
    work->s           = qpalm_calloc(m*2, sizeof(array_element));
    work->index_L     = qpalm_calloc(m*2, sizeof(c_int));
    work->index_P     = qpalm_calloc(m*2, sizeof(c_int));
    work->index_J     = qpalm_calloc(m*2, sizeof(c_int));

    // Primal infeasibility variables
    work->delta_y   = qpalm_calloc(m, sizeof(c_float));
    work->Atdelta_y = qpalm_calloc(n, sizeof(c_float));

    // Dual infeasibility variables
    work->delta_x  = qpalm_calloc(n, sizeof(c_float));
    work->Qdelta_x = qpalm_calloc(n, sizeof(c_float));
    work->Adelta_x = qpalm_calloc(m, sizeof(c_float));

    qpalm_set_factorization_method(work, c);
    // c = &common;
    
    // Allocate scaling structure
    work->scaling       = qpalm_malloc(sizeof(QPALMScaling));
    work->scaling->D    = qpalm_calloc(n, sizeof(c_float));
    work->scaling->Dinv = qpalm_calloc(n, sizeof(c_float));
    work->scaling->E    = qpalm_calloc(m, sizeof(c_float));
    work->scaling->Einv = qpalm_calloc(m, sizeof(c_float));

    work->solver->E_temp = qpalm_calloc(m, sizeof(c_float));
    work->E_temp = work->solver->E_temp;
    work->solver->D_temp = qpalm_calloc(n, sizeof(c_float));
    work->D_temp = work->solver->D_temp;

    // Solver variables
    work->solver->active_constraints = qpalm_calloc(m, sizeof(c_int));
    work->solver->active_constraints_old = qpalm_calloc(m, sizeof(c_int));
    vec_set_scalar_int(work->solver->active_constraints_old, FALSE, m);
    work->solver->reset_newton = TRUE;
    work->solver->enter = qpalm_calloc(m, sizeof(c_int));
    work->solver->leave = qpalm_calloc(m, sizeof(c_int));

    if (work->solver->factorization_method == FACTORIZE_KKT)
    {
        work->solver->rhs_kkt = qpalm_malloc((n+m)*sizeof(c_float));
        work->solver->sol_kkt = qpalm_malloc((n+m)*sizeof(c_float));
        c_int kkt_nzmax = work->data->Q->nzmax + work->data->A->nzmax + m;
        work->solver->kkt_full = ladel_sparse_alloc(n+m, n+m, kkt_nzmax, UPPER, TRUE, FALSE);
        work->solver->kkt = ladel_sparse_alloc(n+m, n+m, kkt_nzmax, UPPER, TRUE, TRUE);
        work->solver->first_row_A = qpalm_malloc(m*sizeof(c_int));
        work->solver->first_elem_A = qpalm_malloc(m*sizeof(c_float));
        work->solver->sym = ladel_symbolics_alloc(m+n);
    } 
    else if (work->solver->factorization_method == FACTORIZE_SCHUR)
    {
        work->solver->sym = ladel_symbolics_alloc(n);
    }
    
    work->solver->neg_dphi = qpalm_calloc(n, sizeof(c_float));
    work->neg_dphi = work->solver->neg_dphi; 
    work->solver->d = qpalm_calloc(n, sizeof(c_float));
    work->d = work->solver->d;
    work->solver->Qd = qpalm_calloc(n, sizeof(c_float));
    work->Qd = work->solver->Qd;
    work->solver->Ad = qpalm_calloc(m, sizeof(c_float));
    work->Ad = work->solver->Ad;
    work->solver->yh = qpalm_calloc(m, sizeof(c_float));
    work->yh = work->solver->yh;
    work->solver->Atyh = qpalm_calloc(n, sizeof(c_float));
    work->Atyh = work->solver->Atyh;
    work->solver->At_scale = qpalm_calloc(m, sizeof(c_float));

    work->solver->first_factorization = TRUE;
    
    if (work->settings->enable_dual_termination)
        work->solver->sym_Q = ladel_symbolics_alloc(n);

    // Allocate solution
    work->solution    = qpalm_calloc(1, sizeof(QPALMSolution));
    work->solution->x = qpalm_calloc(1, n * sizeof(c_float));
    work->solution->y = qpalm_calloc(1, m * sizeof(c_float));
    
    // Allocate and initialize information
    work->info                = qpalm_calloc(1, sizeof(QPALMInfo));
    update_status(work->info, QPALM_UNSOLVED);
    # ifdef QPALM_TIMING
    work->info->solve_time  = 0.0;                    // Solve time to zero
    work->info->run_time    = 0.0;                    // Total run time to zero
    work->info->setup_time  = qpalm_toc(work->timer); // Update timer information
    # endif /* ifdef QPALM_TIMING */
    atomic_init(&work->cancel, 0);

    return work;
}


void qpalm_warm_start(QPALMWorkspace *work, const c_float *x_warm_start, const c_float *y_warm_start) 
{
    atomic_store(&work->cancel, 0);

    // If we have previously solved the problem, then reset the setup time
    if (work->info->status_val != QPALM_UNSOLVED) 
    {
        #ifdef QPALM_TIMING
        work->info->setup_time = 0;
        #endif /* ifdef QPALM_TIMING */
        work->info->status_val = QPALM_UNSOLVED;
    }
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer);
    #endif
    
    size_t n = work->data->n;
    size_t m = work->data->m;
    
    if (x_warm_start != NULL) 
    {
        prea_vec_copy(x_warm_start, work->x, n);
    } 
    else 
    {
        vec_set_scalar(work->x, 0., n);
    }

    if (y_warm_start != NULL) 
    {
        prea_vec_copy(y_warm_start, work->y, m);
    } 
    else 
    {
        vec_set_scalar(work->y, 0., m);
    }
    
    work->initialized = TRUE;

    #ifdef QPALM_TIMING
    work->info->setup_time += qpalm_toc(work->timer); // Start timer
    #endif /* ifdef QPALM_TIMING */

}

static void qpalm_initialize(QPALMWorkspace *work, solver_common **common1, solver_common **common2)
{
    // If we have previously solved the problem, then reset the setup time
    if (work->info->status_val != QPALM_UNSOLVED) 
    {
        #ifdef QPALM_TIMING
        work->info->setup_time = 0;
        #endif /* ifdef QPALM_TIMING */
        work->info->status_val = QPALM_UNSOLVED;
    }
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer);
    #endif

    // Print header
    #ifdef QPALM_PRINTING
    if (work->settings->verbose) 
    {
        print_header();
    }
    #endif

    size_t n = work->data->n;
    size_t m = work->data->m;
    *common1 = ladel_workspace_allocate(n+m);
    if (work->settings->enable_dual_termination) *common2 = ladel_workspace_allocate(n);
    else *common2 = *common1;
    solver_common *c = *common1, *c2 = *common2;

    if (!work->initialized) 
    {
        qpalm_warm_start(work, NULL, NULL);
    } 

    work->eps_abs_in = work->settings->eps_abs_in;
    work->eps_rel_in = work->settings->eps_rel_in;
    work->solver->first_factorization = TRUE;
    work->solver->reset_newton = TRUE;
    work->gamma = work->settings->gamma_init;
    work->gamma_maxed = FALSE;
    vec_set_scalar_int(work->solver->active_constraints_old, FALSE, work->data->m);

    if (work->x == NULL)
    {
        work->x = qpalm_calloc(n, sizeof(c_float));
        vec_set_scalar(work->x, 0., n);
        vec_set_scalar(work->x_prev, 0., n);
        vec_set_scalar(work->x0, 0., n);
        vec_set_scalar(work->Qx, 0., n);
        vec_set_scalar(work->Ax, 0., m);
        work->info->objective = work->data->c;
    }
    else
    {
        mat_vec(work->data->Q, work->x, work->Qx, c);
        mat_vec(work->data->A, work->x, work->Ax, c);
    }

    if (work->y == NULL)
    {
        work->y = qpalm_calloc(m, sizeof(c_float));
        vec_set_scalar(work->y, 0., m);
    }

    for (size_t i = 0; i < work->data->m; i++)
    {
        if (work->data->bmax[i] > QPALM_INFTY) work->data->bmax[i] = QPALM_INFTY;
        if (work->data->bmin[i] < -QPALM_INFTY) work->data->bmin[i] = -QPALM_INFTY;
    }
    
    if (work->settings->scaling) 
    {
        scale_data(work); 
    }

    //Actions to perform after scaling
    prea_vec_copy(work->x, work->x0, n);
    prea_vec_copy(work->x, work->x_prev, n);
    mat_tpose_vec(work->data->A, work->y, work->Aty, c);
    prea_vec_copy(work->y, work->yh, m);
    prea_vec_copy(work->Aty, work->Atyh, n);

    if (work->solver->factorization_method == FACTORIZE_KKT)
    {
        if (work->solver->At) ladel_sparse_free(work->solver->At);
        work->solver->At = ladel_transpose(work->data->A, TRUE, c);
    }

    if (work->settings->nonconvex) 
    {
        set_settings_nonconvex(work, c);
    }
    if (work->settings->proximal) 
    {
        vec_add_scaled(work->Qx, work->x, work->Qx, 1/work->gamma, n);
    } 

    work->info->objective = compute_objective(work);
    
    initialize_sigma(work, c);

    //Provide LD factor of Q in case dual_termination is enabled
    //NB assume Q is positive definite
    if (work->settings->enable_dual_termination) 
    {
        if (work->solver->LD_Q) ladel_factor_free(work->solver->LD_Q);
        ladel_factorize(work->data->Q, work->solver->sym_Q, work->settings->ordering, &work->solver->LD_Q, c2);
        work->info->dual_objective = compute_dual_objective(work, c2);    
    } else 
    {
        work->info->dual_objective = QPALM_NULL;
    }

    #ifdef QPALM_TIMING
    work->info->setup_time += qpalm_toc(work->timer); // Start timer
    #endif /* ifdef QPALM_TIMING */
}

static void qpalm_termination(QPALMWorkspace *work, solver_common* c, solver_common *c2, c_int iter, c_int iter_out)
{
    if (work->info->status_val == QPALM_SOLVED ||
        work->info->status_val == QPALM_DUAL_TERMINATED ||
        work->info->status_val == QPALM_TIME_LIMIT_REACHED ||
        work->info->status_val == QPALM_MAX_ITER_REACHED ||
        work->info->status_val == QPALM_USER_CANCELLATION)
    {
        if (work->info->status_val == QPALM_TIME_LIMIT_REACHED ||
            work->info->status_val == QPALM_MAX_ITER_REACHED ||
            work->info->status_val == QPALM_USER_CANCELLATION)
        {
            compute_residuals(work, c);
        }
        if (work->settings->enable_dual_termination)
        {
            work->info->dual_objective = compute_dual_objective(work, c2);
        }
        store_solution(work);
    }
    else if (work->info->status_val == QPALM_PRIMAL_INFEASIBLE)
    {
        if (work->settings->scaling) 
        {
            vec_self_mult_scalar(work->delta_y, work->scaling->cinv, work->data->m);
            vec_ew_prod(work->scaling->E, work->delta_y, work->delta_y, work->data->m);
        }
    }
    else if (work->info->status_val == QPALM_DUAL_INFEASIBLE)
    {
        if (work->settings->scaling) 
        {
            vec_ew_prod(work->scaling->D, work->delta_x, work->delta_x, work->data->n);
        }
    }

    unscale_data(work);

    work->initialized = FALSE;
    work->info->iter = iter;
    work->info->iter_out = iter_out;

    /* Update solve time and run time */
    #ifdef QPALM_TIMING
        work->info->solve_time = qpalm_toc(work->timer);
        work->info->run_time = work->info->setup_time +
                        work->info->solve_time;
    #endif /* ifdef QPALM_TIMING */
    
    c = ladel_workspace_free(c);
    if (work->settings->enable_dual_termination) 
        c2 = ladel_workspace_free(c2);

    #ifdef QPALM_PRINTING
    if (work->settings->verbose) 
    {
        print_iteration(iter, work); 
        print_final_message(work);
    }
    #endif
}

static void qpalm_terminate_on_status(QPALMWorkspace *work, solver_common *c, solver_common *c2, c_int iter, c_int iter_out, c_int status_val)
{
    update_status(work->info, status_val);
    qpalm_termination(work, c, c2, iter, iter_out);
}

void qpalm_cancel(QPALMWorkspace *work) {
    atomic_store(&work->cancel, 1);
}

void qpalm_solve(QPALMWorkspace *work) 
{
	#if defined(QPALM_PRINTING) && defined(_WIN32) && defined(_MSC_VER) && _MSC_VER < 1900
	unsigned int print_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
	#endif

    //Initialize ladel workspace, qpalm variables and perform scaling (timing added to setup)
    solver_common *c, *c2;
    qpalm_initialize(work, &c, &c2);

    // Start the timer for the solve routine
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer); 
    c_float current_time;
    #endif /* ifdef QPALM_TIMING */
    
    size_t m = work->data->m;
    c_int iter;
    c_int iter_out = 0;
    c_int prev_iter = 0; /* iteration number at which the previous subproblem finished*/
    c_float eps_k_abs = work->settings->eps_abs_in; 
    c_float eps_k_rel = work->settings->eps_rel_in; 
    c_int no_change_in_active_constraints = 0;

    for (iter = 0; iter < work->settings->max_iter; iter++) 
    {
        /* Check whether we passed the time limit */
        #ifdef QPALM_TIMING
        current_time = work->info->setup_time + qpalm_toc(work->timer); // Start timer
        if (current_time > work->settings->time_limit) 
        {
            qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_TIME_LIMIT_REACHED);
            return;
        }
        #endif /* ifdef QPALM_TIMING */
        if (atomic_load(&work->cancel))
        {
            qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_USER_CANCELLATION);
            return;
        }

        /*Perform the iteration */
        compute_residuals(work, c);
        calculate_residual_norms_and_tolerances(work);
        
        if (is_solved(work)) 
        {
            qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_SOLVED);
            return;
        }
        else if (is_primal_infeasible(work)) 
        {
            qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_PRIMAL_INFEASIBLE);
            return;
        } 
        else if (is_dual_infeasible(work)) 
        {
            qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_DUAL_INFEASIBLE);
            return;
        }
        else if (check_subproblem_termination(work) || (no_change_in_active_constraints == 3)) 
        {
            update_dual_iterate_and_parameters(work, c, iter_out, &eps_k_abs, &eps_k_rel);
        
            if(work->settings->enable_dual_termination) 
            {
                work->info->dual_objective = compute_dual_objective(work, c);
                if (work->info->dual_objective > work->settings->dual_objective_limit) 
                {
                    qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_DUAL_TERMINATED);
                    return; 
                }
            }

            no_change_in_active_constraints = 0;
            iter_out++;
            prev_iter = iter;

            #ifdef QPALM_PRINTING
            if (work->settings->verbose && mod(iter, work->settings->print_iter) == 0) 
            {
                qpalm_print("%4" LADEL_PRIi " | ---------------------------------------------------\n", iter);
            }
            #endif    
        } 
        else if (iter == prev_iter + work->settings->inner_max_iter) /*subproblem is hanging so try updating params*/
        { 
            if (iter_out > 0 && work->info->pri_res_norm > work->eps_pri) 
            {
                update_sigma(work, c);
            } 

            if(work->settings->proximal) 
            {
                update_gamma(work);
                if (!work->settings->nonconvex) prea_vec_copy(work->x, work->x0, work->data->n);
            }

            prea_vec_copy(work->pri_res, work->pri_res_in, m);

            no_change_in_active_constraints = 0;     
            iter_out++;
            prev_iter = iter;
        } 
        else /*primal update*/
        {
            if (work->solver->nb_enter + work->solver->nb_leave) no_change_in_active_constraints = 0;
            else no_change_in_active_constraints++;

            if (mod(iter, work->settings->reset_newton_iter) == 0) work->solver->reset_newton = TRUE; 
            update_primal_iterate(work, c);

            #ifdef QPALM_PRINTING
            if (work->settings->verbose && mod(iter, work->settings->print_iter) == 0) 
            {
                work->info->objective = compute_objective(work);
                print_iteration(iter, work);
            }
            #endif
        }
    }

    /* If we get here, qpalm has unfortunately hit the maximum number of iterations */
    qpalm_terminate_on_status(work, c, c2, iter, iter_out, QPALM_MAX_ITER_REACHED);
    return;
}

void qpalm_update_settings(QPALMWorkspace* work, const QPALMSettings *settings) 
{
    // If we have previously solved the problem, then reset the setup time
    if (work->info->status_val != QPALM_UNSOLVED) 
    {
        #ifdef QPALM_TIMING
        work->info->setup_time = 0;
        #endif /* ifdef QPALM_TIMING */
        work->info->status_val = QPALM_UNSOLVED;
    }
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer); // Start timer
    #endif /* ifdef QPALM_TIMING */
    
    // Validate settings
    if (!validate_settings(settings)) 
    {
        # ifdef QPALM_PRINTING
        qpalm_eprint("Settings validation returned failure");
        # endif /* ifdef QPALM_PRINTING */
        update_status(work->info, QPALM_ERROR);
        return;
    }
    
    // Copy settings
    qpalm_free(work->settings);
    work->settings = copy_settings(settings);
    work->sqrt_delta = c_sqrt(work->settings->delta);
    # ifdef QPALM_TIMING
    work->info->setup_time += qpalm_toc(work->timer);
    # endif /* ifdef QPALM_TIMING */
}

void qpalm_update_bounds(QPALMWorkspace *work, const c_float *bmin, const c_float *bmax) 
{
    // If we have previously solved the problem, then reset the setup time
    if (work->info->status_val != QPALM_UNSOLVED) 
    {
        #ifdef QPALM_TIMING
        work->info->setup_time = 0;
        #endif /* ifdef QPALM_TIMING */
        work->info->status_val = QPALM_UNSOLVED;
    }
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer); // Start timer
    #endif /* ifdef QPALM_TIMING */

    // Validate bounds
    size_t j;
    size_t m = work->data->m;
    if (bmin != NULL && bmax != NULL) 
    {
        for (j = 0; j < m; j++) 
        {
            if (bmin[j] > bmax[j]) 
            {
                # ifdef QPALM_PRINTING
                qpalm_eprint("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
                        (int)j, work->data->bmin[j], work->data->bmax[j]);
                # endif /* ifdef QPALM_PRINTING */
                update_status(work->info, QPALM_ERROR);
                return;
            }
        }
    }
    
    if (bmin != NULL) 
    {
        prea_vec_copy(bmin, work->data->bmin, m); 
    }
    if (bmax != NULL) 
    {
        prea_vec_copy(bmax, work->data->bmax, m);
    }     
    
    # ifdef QPALM_TIMING
    work->info->setup_time += qpalm_toc(work->timer);
    # endif /* ifdef QPALM_TIMING */
}

void qpalm_update_q(QPALMWorkspace *work, const c_float *q) 
{
    // If we have previously solved the problem, then reset the setup time
    if (work->info->status_val != QPALM_UNSOLVED) 
    {
        #ifdef QPALM_TIMING
        work->info->setup_time = 0;
        #endif /* ifdef QPALM_TIMING */
        work->info->status_val = QPALM_UNSOLVED;
    }
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer); // Start timer
    #endif /* ifdef QPALM_TIMING */

    size_t n = work->data->n;
    prea_vec_copy(q, work->data->q, n);    
    # ifdef QPALM_TIMING
    work->info->setup_time += qpalm_toc(work->timer);
    # endif /* ifdef QPALM_TIMING */
}

void qpalm_update_Q_A(QPALMWorkspace *work, const c_float *Qx, const c_float *Ax)
{
    work->solver->first_factorization = TRUE;
    // If we have previously solved the problem, then reset the setup time
    if (work->info->status_val != QPALM_UNSOLVED) 
    {
        #ifdef QPALM_TIMING
        work->info->setup_time = 0;
        #endif /* ifdef QPALM_TIMING */
        work->info->status_val = QPALM_UNSOLVED;
    }
    #ifdef QPALM_TIMING
    qpalm_tic(work->timer); // Start timer
    #endif /* ifdef QPALM_TIMING */

    ladel_sparse_matrix *Q = work->data->Q, *A = work->data->A;
    prea_vec_copy(Qx, Q->x, Q->nzmax);
    prea_vec_copy(Ax, A->x, A->nzmax);
    
    # ifdef QPALM_TIMING
    work->info->setup_time += qpalm_toc(work->timer);
    # endif /* ifdef QPALM_TIMING */
}

void qpalm_cleanup(QPALMWorkspace *work) 
{
    if (work) 
    { 
        // Free Data
        if (work->data) 
        {
            work->data->Q = ladel_sparse_free(work->data->Q);

            work->data->A = ladel_sparse_free(work->data->A);

            if (work->data->q) qpalm_free(work->data->q);

            if (work->data->bmin) qpalm_free(work->data->bmin);

            if (work->data->bmax) qpalm_free(work->data->bmax);
            qpalm_free(work->data);
        }

        // Free scaling
        if (work->scaling->D) qpalm_free(work->scaling->D);

        if (work->scaling->Dinv) qpalm_free(work->scaling->Dinv);

        if (work->scaling->E) qpalm_free(work->scaling->E);

        if (work->scaling->Einv) qpalm_free(work->scaling->Einv);

        qpalm_free(work->scaling);

        // Free other Variables
        if (work->x) qpalm_free(work->x);
        
        if (work->y) qpalm_free(work->y);

        if (work->Ax) qpalm_free(work->Ax);

        if (work->Qx) qpalm_free(work->Qx);
        
        if (work->x_prev) qpalm_free(work->x_prev);

        if (work->Aty) qpalm_free(work->Aty);

        if (work->temp_m) qpalm_free(work->temp_m);

        if (work->temp_n) qpalm_free(work->temp_n);

        if (work->sigma) qpalm_free(work->sigma);
        
        if (work->sigma_inv) qpalm_free(work->sigma_inv);

        if (work->z) qpalm_free(work->z);

        if (work->Axys) qpalm_free(work->Axys);

        if (work->pri_res) qpalm_free(work->pri_res);

        if (work->pri_res_in) qpalm_free(work->pri_res_in);

        if (work->df) qpalm_free(work->df);

        if (work->x0) qpalm_free(work->x0);

        if (work->xx0) qpalm_free(work->xx0);

        if (work->dphi) qpalm_free(work->dphi);

        if (work->dphi_prev) qpalm_free(work->dphi_prev);

        if (work->sqrt_sigma) qpalm_free(work->sqrt_sigma);

        if (work->delta) qpalm_free(work->delta);

        if (work->alpha) qpalm_free(work->alpha);

        if (work->delta2) qpalm_free(work->delta2);

        if (work->delta_alpha) qpalm_free(work->delta_alpha);

        if (work->temp_2m) qpalm_free(work->temp_2m);

        if (work->s) qpalm_free(work->s);

        if (work->index_L) qpalm_free(work->index_L);

        if (work->index_P) qpalm_free(work->index_P);

        if (work->index_J) qpalm_free(work->index_J);

        if (work->delta_y) qpalm_free(work->delta_y);

        if (work->Atdelta_y) qpalm_free(work->Atdelta_y);

        if (work->delta_x) qpalm_free(work->delta_x);

        if (work->Qdelta_x) qpalm_free(work->Qdelta_x);

        if (work->Adelta_x) qpalm_free(work->Adelta_x);

        // Free Settings
        if (work->settings) qpalm_free(work->settings);

        //Free chol struct
        if (work->solver) 
        {
            if (work->solver->active_constraints) qpalm_free(work->solver->active_constraints);

            if (work->solver->active_constraints_old) qpalm_free(work->solver->active_constraints_old);

            if (work->solver->enter) qpalm_free(work->solver->enter);

            if (work->solver->leave) qpalm_free(work->solver->leave);


            work->solver->sol_kkt = ladel_free(work->solver->sol_kkt);

            work->solver->rhs_kkt = ladel_free(work->solver->rhs_kkt);

            work->solver->D_temp = ladel_free(work->solver->D_temp);

            work->solver->E_temp = ladel_free(work->solver->E_temp);

            work->solver->neg_dphi = ladel_free(work->solver->neg_dphi);

            work->solver->d = ladel_free(work->solver->d);

            work->solver->Qd = ladel_free(work->solver->Qd);

            work->solver->Ad = ladel_free(work->solver->Ad);

            work->solver->yh = ladel_free(work->solver->yh);

            work->solver->Atyh = ladel_free(work->solver->Atyh);

            work->solver->LD = ladel_factor_free(work->solver->LD);

            work->solver->LD_Q = ladel_factor_free(work->solver->LD_Q);

            work->solver->sym = ladel_symbolics_free(work->solver->sym);

            work->solver->sym_Q = ladel_symbolics_free(work->solver->sym_Q);

            work->solver->At_scale = ladel_free(work->solver->At_scale);

            work->solver->At_sqrt_sigma = ladel_sparse_free(work->solver->At_sqrt_sigma);

            work->solver->At = ladel_sparse_free(work->solver->At);

            work->solver->kkt = ladel_sparse_free(work->solver->kkt);

            work->solver->kkt_full = ladel_sparse_free(work->solver->kkt_full);

            work->solver->first_row_A = ladel_free(work->solver->first_row_A);

            work->solver->first_elem_A = ladel_free(work->solver->first_elem_A);


            qpalm_free(work->solver);      
        }
        
        // Free solution
        if (work->solution) 
        {
            if (work->solution->x) qpalm_free(work->solution->x);

            if (work->solution->y) qpalm_free(work->solution->y);
            qpalm_free(work->solution);
        }

        // Free timer
        # ifdef QPALM_TIMING
        if (work->timer) qpalm_free(work->timer);
        # endif /* ifdef QPALM_TIMING */

        // Free information
        if (work->info) qpalm_free(work->info);

        // Free work
        qpalm_free(work);
    }

}


# ifdef __cplusplus
}
# endif // ifdef __cplusplus
