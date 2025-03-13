#include <mex.h>
#include <string.h>
#include <qpalm/global_opts.h>
#include <qpalm.h>
#include <qpalm/util.h>
#include <qpalm/constants.h>
#include <ladel.h>
#include <ladel_mex_util.h>

//Modes of operation
#define MODE_DEFAULT_SETTINGS "default_settings"
#define MODE_RETRIEVE_CONSTANT "constant"
#define MODE_SETUP "setup"
#define MODE_WARM_START "warm_start"
#define MODE_UPDATE_BOUNDS "update_bounds"
#define MODE_UPDATE_LINEAR "update_q"
#define MODE_UPDATE_MATRICES "update_Q_A"
#define MODE_SOLVE "solve"
#define MODE_DELETE "delete"

// QPALM work identifier
static QPALMWorkspace* qpalm_work = NULL;  

// all of the QPALM_INFO fieldnames as strings
const char* QPALM_INFO_FIELDS[] = {"iter",          //c_int
                                  "iter_out",       //c_int
                                  "status" ,        //char*
                                  "status_val" ,    //c_int
                                  "pri_res_norm",   //c_float
                                  "dua_res_norm",   //c_float
                                  "dua2_res_norm",  //c_float
                                  "objective",      //c_float
                                  "dual_objective", //c_float
                                  "setup_time",     //c_float, only used if QPALM_TIMING
                                  "solve_time",     //c_float, only used if QPALM_TIMING
                                  "run_time"};      //c_float, only used if QPALM_TIMING


const char* QPALM_SETTINGS_FIELDS[] = {"max_iter",                  //c_int
                                      "inner_max_iter",             //c_int  
                                      "eps_abs",                    //c_float
                                      "eps_rel",                    //c_float
                                      "eps_abs_in",                 //c_float
                                      "eps_rel_in",                 //c_float
                                      "rho",                        //c_float
                                      "eps_prim_inf",               //c_float
                                      "eps_dual_inf",               //c_float
                                      "theta",                      //c_float
                                      "delta",                      //c_float
                                      "sigma_max",                  //c_float
                                      "sigma_init",                 //c_float
                                      "proximal",                   //c_int
                                      "gamma_init",                 //c_float
                                      "gamma_upd",                  //c_float
                                      "gamma_max",                  //c_float
                                      "scaling",                    //c_int
                                      "nonconvex",                  //c_int
                                      "warm_start",                 //c_int
                                      "print_iter",                 //c_int
                                      "reset_newton_iter",          //c_int
                                      "enable_dual_termination",    //c_int
                                      "dual_objective_limit",       //c_int
                                      "time_limit",                 //c_float
                                      "verbose",                    //c_int
                                      "ordering",                   //c_int
                                      "factorization_method",       //c_int    
                                      "max_rank_update",            //c_int
                                      "max_rank_update_fraction"};  //c_float                      


// internal utility functions
void      castToDoubleArr(c_float *arr, double* arr_out, size_t len);
void      setToNaN(double* arr_out, size_t len);
void      copyMxStructToSettings(const mxArray*, QPALMSettings*);
mxArray*  copySettingsToMxStruct(QPALMSettings* settings);
mxArray*  copyInfoToMxStruct(QPALMInfo* info);

/* Use Mex's custom allocators */
static void *qpalm_mex_calloc(size_t n, size_t size) {
    void *m = mxCalloc(n, size);
    mexMakeMemoryPersistent(m);
    return m;
}
static void *qpalm_mex_malloc(size_t size) {
    void *m = mxMalloc(size);
    mexMakeMemoryPersistent(m);
    return m;
}
static void* qpalm_mex_realloc(void *p, size_t size) {
    void *p_new = mxRealloc(p, size);
    mexMakeMemoryPersistent(p_new);
    return p_new;
}
static void qpalm_mex_free(void* p) {
    mxFree(p);
}

struct {
    calloc_sig *calloc;
    malloc_sig *malloc;
    realloc_sig *realloc;
    free_sig *free;
    printf_sig *printf;
} static old_config = {
    .calloc = NULL,
    .malloc = NULL,
    .realloc = NULL,
    .free = NULL,
    .printf = NULL
};

static int config_configured = 0;

static void configure_mex_specific_functions(void) {
    if (config_configured)
        return;
    config_configured = 1;
    old_config.calloc = qpalm_set_alloc_config_calloc(&qpalm_mex_calloc);
    old_config.malloc = qpalm_set_alloc_config_malloc(&qpalm_mex_malloc);
    old_config.realloc = qpalm_set_alloc_config_realloc(&qpalm_mex_realloc);
    old_config.free = qpalm_set_alloc_config_free(&qpalm_mex_free);
    old_config.printf = qpalm_set_print_config_printf(&mexPrintf);
}

static void unconfigure_mex_specific_functions(void) {
    if (old_config.calloc) {
        qpalm_set_alloc_config_calloc(old_config.calloc);
        old_config.calloc = NULL;
    }
    if (old_config.malloc) {
        qpalm_set_alloc_config_malloc(old_config.malloc);
        old_config.malloc = NULL;
    }
    if (old_config.realloc) {
        qpalm_set_alloc_config_realloc(old_config.realloc);
        old_config.realloc = NULL;
    }
    if (old_config.free) {
        qpalm_set_alloc_config_free(old_config.free);
        old_config.free = NULL;
    }
    if (old_config.printf) {
        qpalm_set_print_config_printf(old_config.printf);
        old_config.printf = NULL;
    }
}

/**
 * Function that mex calls when it closes unexpectedly
 * Frees the workspace
 */

void exitFcn() {
  if (qpalm_work != NULL) {
      qpalm_cleanup(qpalm_work);
      qpalm_work = NULL;
  }
  if (config_configured) {
      unconfigure_mex_specific_functions();
  }
}

/**
 * The gateway function to QPALM
 *
 * qpalm_mex('default_settings');
 * qpalm_mex('constant', constant_name)
 * qpalm_mex('setup',n,m,Q,q,A,bmin,bmax,theSettings);
 * qpalm_mex('warm_start', x, y);
 * [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = qpalm_mex('solve');
 * qpalm_mex('delete');
 *
 * @param nlhs Number of output arguments
 * @param plhs Array of output argument pointers
 * @param nrhs Number of input arguments
 * @param prhs Array of input argument pointers
 */
void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {
    
    // Set function to call when mex closes unexpectedly
    mexAtExit(exitFcn);

    /* Set the custom allocator and print functions */
    configure_mex_specific_functions();

    // Get the command string
    char cmd[64];

    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    
    // report the default settings
    if (strcmp(cmd, MODE_DEFAULT_SETTINGS) == 0) {
        // Warn if other commands were ignored
        if (nrhs > 2)
            mexWarnMsgTxt("Default settings: unexpected number of arguments.");

        //Create a Settings structure in default form and report the results
        //Useful for external solver packages (e.g. Yalmip) that want to
        //know which solver settings are supported
        QPALMSettings* defaults = (QPALMSettings *)mxCalloc(1,sizeof(QPALMSettings));
        qpalm_set_default_settings(defaults);
        plhs[0] = copySettingsToMxStruct(defaults);
        mxFree(defaults);
        return;

    } else if (strcmp(cmd, MODE_DELETE) == 0) {
        
        //clean up the problem workspace
        if(qpalm_work != NULL){
            qpalm_cleanup(qpalm_work);
            qpalm_work = NULL;
        }
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 1)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;

    } else if (strcmp(cmd, MODE_SETUP) == 0) {

        //throw an error if this is called more than once
        if(qpalm_work != NULL){
            mexErrMsgTxt("Solver is already initialized with problem data.");
        }   

        //Create data and settings containers
        QPALMSettings* settings = (QPALMSettings *)mxCalloc(1,sizeof(QPALMSettings));
        QPALMData*     data     = (QPALMData *)mxCalloc(1,sizeof(QPALMData));

        // handle the problem data first.  Matlab-side
        // class wrapper is responsible for ensuring that
        // Q and A are sparse matrices,  everything
        // else is a dense vector and all inputs are
        // of compatible dimension

        // Get mxArrays
        const mxArray* Q  = prhs[3];
        const mxArray* q  = prhs[4];
        const mxArray* A  = prhs[5];
        const mxArray* bmin = prhs[6];
        const mxArray* bmax = prhs[7];

        // Create Data Structure
        data->n = (size_t) mxGetScalar(prhs[1]);
        data->m = (size_t) mxGetScalar(prhs[2]);
        data->c = 0; //TODO: Allow for a constant to be passed.
        data->q = mxGetPr(q);
        data->bmin = mxGetPr(bmin);
        data->bmax = mxGetPr(bmax);

        // Convert matrices from matlab to C structure (LADEL/CHOLMOD)
        solver_sparse Amatrix, Qmatrix;
        data->A = ladel_get_sparse_from_matlab(A, &Amatrix, UNSYMMETRIC);
        data->Q = ladel_get_sparse_from_matlab(Q, &Qmatrix, UPPER);

        // Create Settings
        const mxArray* mxSettings = prhs[8];
        if(mxIsEmpty(mxSettings)){
            // use defaults
            qpalm_set_default_settings(settings);
        } else {
            //populate settings structure from mxArray input
            copyMxStructToSettings(mxSettings, settings);
        }
        
        // Setup workspace
        qpalm_work = qpalm_setup(data, settings);

        if(qpalm_work == NULL){
            mexErrMsgTxt("Invalid problem setup");
        }

         // cleanup temporary structures
         // Don't free data->q, data->bmin, data->bmin because they are pointers to the mxArrays
         // Don't free data->A and data->Q because they are only a shallow copy
         mxFree(data);
         mxFree(settings);

        return;


    } else if (strcmp(cmd, MODE_WARM_START) == 0) {
        if (nlhs != 0 || nrhs != 3 ){
            mexErrMsgTxt("Solve : wrong number of inputs / outputs");
        }
        if(!qpalm_work){
            mexErrMsgTxt("Work is not setup.");
        }

        c_float *x, *y;

        const mxArray* xmatlab  = prhs[1];
        const mxArray* ymatlab  = prhs[2];

        if (mxIsEmpty(xmatlab)) {
            x = NULL;
        } else {
            x = mxGetPr(xmatlab);
        }
            
        if (mxIsEmpty(ymatlab)) {
            y = NULL;
        } else {
            y = mxGetPr(ymatlab);
        }

        qpalm_warm_start(qpalm_work, x, y);

        return;

    } else if (strcmp(cmd, MODE_UPDATE_BOUNDS) == 0) {
        
        if (nlhs != 0 || nrhs != 3){
            mexErrMsgTxt("Update bounds : wrong number of inputs / outputs");
        }
        if(!qpalm_work){
            mexErrMsgTxt("Work is not setup.");
        }
        
        c_float *bmin, *bmax;

        const mxArray* bmin_matlab  = prhs[1];
        const mxArray* bmax_matlab  = prhs[2];

        if (mxIsEmpty(bmin_matlab)) {
            bmin = NULL;
        } else {
            bmin = mxGetPr(bmin_matlab);
        }
            
        if (mxIsEmpty(bmax_matlab)) {
            bmax = NULL;
        } else {
            bmax = mxGetPr(bmax_matlab);
        }

        qpalm_update_bounds(qpalm_work, bmin, bmax);

    } else if (strcmp(cmd, MODE_UPDATE_LINEAR) == 0) {
        
        if (nlhs != 0 || nrhs != 2){
            mexErrMsgTxt("Update q : wrong number of inputs / outputs");
        }
        if(!qpalm_work){
            mexErrMsgTxt("Work is not setup.");
        }
        
        if (!mxIsEmpty(prhs[1])) {
            c_float *q = mxGetPr(prhs[1]);
            qpalm_update_q(qpalm_work, q);
        } else {
            mexWarnMsgTxt("Update q: Empty q has no effect.");
        }

    } else if (strcmp(cmd, MODE_UPDATE_MATRICES) == 0) {
        
        if (nlhs != 0 || nrhs != 3){
            mexErrMsgTxt("Update Q and A : wrong number of inputs / outputs");
        }
        if(!qpalm_work){
            mexErrMsgTxt("Work is not setup.");
        }
        
        if (!mxIsEmpty(prhs[1])) {
            qpalm_update_Q_A(qpalm_work, mxGetPr(prhs[1]), mxGetPr(prhs[2]));
        } else {
            mexWarnMsgTxt("Update Q and A: Empty Q has no effect.");
        }

    } else if (strcmp(cmd, MODE_SOLVE) == 0) { // SOLVE

        if (nlhs != 5 || nrhs != 1){
            mexErrMsgTxt("Solve : wrong number of inputs / outputs");
        }
        if(!qpalm_work){
            mexErrMsgTxt("Work is not setup.");
        }
        
        qpalm_solve(qpalm_work);

        // Allocate space for solution
        // primal variables
        plhs[0] = mxCreateDoubleMatrix(qpalm_work->data->n,1,mxREAL);
        // dual variables
        plhs[1] = mxCreateDoubleMatrix(qpalm_work->data->m,1,mxREAL);
        // primal infeasibility certificate
        plhs[2] = mxCreateDoubleMatrix(qpalm_work->data->m,1,mxREAL);
        // dual infeasibility certificate
        plhs[3] = mxCreateDoubleMatrix(qpalm_work->data->n,1,mxREAL);

        //copy results to mxArray outputs
        //assume that five outputs will always
        //be returned to matlab-side class wrapper
        if ((qpalm_work->info->status_val != QPALM_PRIMAL_INFEASIBLE) &&
            (qpalm_work->info->status_val != QPALM_DUAL_INFEASIBLE)){

            //primal and dual solutions
            castToDoubleArr(qpalm_work->solution->x, mxGetPr(plhs[0]), qpalm_work->data->n);
            castToDoubleArr(qpalm_work->solution->y, mxGetPr(plhs[1]), qpalm_work->data->m);

            //infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), qpalm_work->data->m);
            setToNaN(mxGetPr(plhs[3]), qpalm_work->data->n);

        } else if (qpalm_work->info->status_val == QPALM_PRIMAL_INFEASIBLE){ //primal infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpalm_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpalm_work->data->m);

            //primal infeasibility certificates
            castToDoubleArr(qpalm_work->delta_y, mxGetPr(plhs[2]), qpalm_work->data->m);

            //dual infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[3]), qpalm_work->data->n);

        } else if (qpalm_work->info->status_val == QPALM_DUAL_INFEASIBLE) { //dual infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpalm_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpalm_work->data->m);

            //primal infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), qpalm_work->data->m);

            //dual infeasibility certificates
            castToDoubleArr(qpalm_work->delta_x, mxGetPr(plhs[3]), qpalm_work->data->n);

        }

        plhs[4] = copyInfoToMxStruct(qpalm_work->info); // Info structure
        return;



    } else if (strcmp(cmd, MODE_RETRIEVE_CONSTANT) == 0) { // Return solver constants
         
        char constant[32];
        mxGetString(prhs[1], constant, sizeof(constant));

        if (!strcmp("QPALM_INFTY", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_INFTY);
            return;
        }
        if (!strcmp("QPALM_NAN", constant)){
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
            return;
        }

        if (!strcmp("QPALM_SOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_SOLVED);
            return;
        }

        if (!strcmp("QPALM_UNSOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_UNSOLVED);
            return;
        }

        if (!strcmp("QPALM_PRIMAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_PRIMAL_INFEASIBLE);
            return;
        }

        if (!strcmp("QPALM_DUAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_DUAL_INFEASIBLE);
            return;
        }

        if (!strcmp("QPALM_MAX_ITER_REACHED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_MAX_ITER_REACHED);
            return;
        }

        mexErrMsgTxt("Constant not recognized.");

        return;


    } else {
        mexErrMsgTxt("Invalid QPALM mode");
    }

}

void castToDoubleArr(c_float *arr, double* arr_out, size_t len) {
    for (size_t i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
}

void setToNaN(double* arr_out, size_t len){
    size_t i;
    for (i = 0; i < len; i++) {
        arr_out[i] = mxGetNaN();
    }
}

mxArray* copyInfoToMxStruct(QPALMInfo* info){

  //create mxArray with the right number of fields
  int nfields  = sizeof(QPALM_INFO_FIELDS) / sizeof(QPALM_INFO_FIELDS[0]);
  mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,QPALM_INFO_FIELDS);

  //map the QPALM_INFO fields one at a time into mxArrays
  //matlab all numeric values as doubles
  mxSetField(mxPtr, 0, "iter",              mxCreateDoubleScalar(info->iter));
  mxSetField(mxPtr, 0, "iter_out",          mxCreateDoubleScalar(info->iter_out));
  mxSetField(mxPtr, 0, "status",            mxCreateString(info->status));
  mxSetField(mxPtr, 0, "status_val",        mxCreateDoubleScalar(info->status_val));
  mxSetField(mxPtr, 0, "pri_res_norm",      mxCreateDoubleScalar(info->pri_res_norm));
  mxSetField(mxPtr, 0, "dua_res_norm",      mxCreateDoubleScalar(info->dua_res_norm));
  mxSetField(mxPtr, 0, "dua2_res_norm",     mxCreateDoubleScalar(info->dua2_res_norm));
  mxSetField(mxPtr, 0, "objective",         mxCreateDoubleScalar(info->objective));
  mxSetField(mxPtr, 0, "dual_objective",    mxCreateDoubleScalar(info->dual_objective));

  #ifdef QPALM_TIMING
  //if not profiling, these fields will be empty
  mxSetField(mxPtr, 0, "setup_time",  mxCreateDoubleScalar(info->setup_time));
  mxSetField(mxPtr, 0, "solve_time",  mxCreateDoubleScalar(info->solve_time));
  mxSetField(mxPtr, 0, "run_time",    mxCreateDoubleScalar(info->run_time));
  #endif

  return mxPtr;

}

mxArray* copySettingsToMxStruct(QPALMSettings* settings){

  int nfields  = sizeof(QPALM_SETTINGS_FIELDS) / sizeof(QPALM_SETTINGS_FIELDS[0]);
  mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,QPALM_SETTINGS_FIELDS);

  //map the QPALM_SETTINGS fields one at a time into mxArrays
  //matlab handles everything as a double
  mxSetField(mxPtr, 0, "max_iter",                  mxCreateDoubleScalar(settings->max_iter));
  mxSetField(mxPtr, 0, "inner_max_iter",            mxCreateDoubleScalar(settings->inner_max_iter));
  mxSetField(mxPtr, 0, "eps_abs",                   mxCreateDoubleScalar(settings->eps_abs));
  mxSetField(mxPtr, 0, "eps_rel",                   mxCreateDoubleScalar(settings->eps_rel));
  mxSetField(mxPtr, 0, "eps_abs_in",                mxCreateDoubleScalar(settings->eps_abs_in));
  mxSetField(mxPtr, 0, "eps_rel_in",                mxCreateDoubleScalar(settings->eps_rel_in));
  mxSetField(mxPtr, 0, "rho",                       mxCreateDoubleScalar(settings->rho));
  mxSetField(mxPtr, 0, "eps_prim_inf",              mxCreateDoubleScalar(settings->eps_prim_inf));
  mxSetField(mxPtr, 0, "eps_dual_inf",              mxCreateDoubleScalar(settings->eps_dual_inf));
  mxSetField(mxPtr, 0, "theta",                     mxCreateDoubleScalar(settings->theta));
  mxSetField(mxPtr, 0, "delta",                     mxCreateDoubleScalar(settings->delta));
  mxSetField(mxPtr, 0, "sigma_max",                 mxCreateDoubleScalar(settings->sigma_max));
  mxSetField(mxPtr, 0, "sigma_init",                mxCreateDoubleScalar(settings->sigma_init));
  mxSetField(mxPtr, 0, "proximal",                  mxCreateDoubleScalar(settings->proximal));
  mxSetField(mxPtr, 0, "gamma_init",                mxCreateDoubleScalar(settings->gamma_init));
  mxSetField(mxPtr, 0, "gamma_upd",                 mxCreateDoubleScalar(settings->gamma_upd));
  mxSetField(mxPtr, 0, "gamma_max",                 mxCreateDoubleScalar(settings->gamma_max));
  mxSetField(mxPtr, 0, "scaling",                   mxCreateDoubleScalar(settings->scaling));
  mxSetField(mxPtr, 0, "nonconvex",                 mxCreateDoubleScalar(settings->nonconvex));
  mxSetField(mxPtr, 0, "warm_start",                mxCreateDoubleScalar(settings->warm_start));
  mxSetField(mxPtr, 0, "verbose",                   mxCreateDoubleScalar(settings->verbose));
  mxSetField(mxPtr, 0, "print_iter",                mxCreateDoubleScalar(settings->print_iter));
  mxSetField(mxPtr, 0, "reset_newton_iter",         mxCreateDoubleScalar(settings->reset_newton_iter));
  mxSetField(mxPtr, 0, "enable_dual_termination",   mxCreateDoubleScalar(settings->enable_dual_termination));
  mxSetField(mxPtr, 0, "dual_objective_limit",      mxCreateDoubleScalar(settings->dual_objective_limit));
  mxSetField(mxPtr, 0, "time_limit",                mxCreateDoubleScalar(settings->time_limit));
  mxSetField(mxPtr, 0, "ordering",                  mxCreateDoubleScalar(settings->ordering));
  mxSetField(mxPtr, 0, "factorization_method",      mxCreateDoubleScalar(settings->factorization_method));
  mxSetField(mxPtr, 0, "max_rank_update",           mxCreateDoubleScalar(settings->max_rank_update));
  mxSetField(mxPtr, 0, "max_rank_update_fraction",  mxCreateDoubleScalar(settings->max_rank_update_fraction));

  return mxPtr;
}


void copyMxStructToSettings(const mxArray* mxPtr, QPALMSettings* settings){

  //this function assumes that only a complete and validated structure
  //will be passed.  matlab mex-side function is responsible for checking
  //structure validity

  //map the QPALM_SETTINGS fields one at a time into mxArrays
  //matlab handles everything as a double
  settings->max_iter                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
  settings->inner_max_iter            = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "inner_max_iter"));
  settings->eps_abs                   = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
  settings->eps_rel                   = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel"));
  settings->eps_abs_in                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs_in"));
  settings->eps_rel_in                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel_in"));
  settings->rho                       = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "rho"));
  settings->eps_prim_inf              = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_prim_inf"));
  settings->eps_dual_inf              = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  settings->theta                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "theta"));
  settings->delta                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->sigma_max                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "sigma_max"));
  settings->sigma_init                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "sigma_init"));
  settings->proximal                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "proximal"));
  settings->gamma_init                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "gamma_init"));  
  settings->gamma_upd                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "gamma_upd"));
  settings->gamma_max                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "gamma_max"));
  settings->scaling                   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "scaling"));
  settings->nonconvex                 = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "nonconvex"));
  settings->warm_start                = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "warm_start"));
  settings->verbose                   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
  settings->print_iter                = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "print_iter"));
  settings->reset_newton_iter         = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "reset_newton_iter"));
  settings->enable_dual_termination   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "enable_dual_termination"));
  settings->dual_objective_limit      = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "dual_objective_limit"));
  settings->time_limit                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "time_limit"));
  settings->ordering                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "ordering"));
  settings->factorization_method      = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "factorization_method"));
  settings->max_rank_update           = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "max_rank_update"));
  settings->max_rank_update_fraction  = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "max_rank_update_fraction"));
}
