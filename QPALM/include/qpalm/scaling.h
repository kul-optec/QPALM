/**
 * @file scaling.h
 * @author Ben Hermans
 * @brief Problem data scaling during setup.
 * @details This file includes the routine that is called during setup to scale the problem data, 
 * and initial guesses if the problem is warm-started. Scaling the problem is useful to prevent 
 * large changes in the active set and to guard against ill-conditioning in the objective function.
 *  
 */

#ifndef SCALING_H
# define SCALING_H

# ifdef __cplusplus
extern "C" {
# endif 

// Functions to scale problem data
# include "types.h"
# include "lin_alg.h"
# include "constants.h"


/**
 * Scale problem matrices
 * 
 * Ruiz scaling \cite ruiz2001scaling is applied to the constraint matrix A. This means that the rows and columns of A are 
 * scaled elementwise by the square root of their infinity norm, and this for a number of work->settings->scaling
 * iterations. The resulting scaling can be written as @f$\bar{A}\leftarrow EAD@f$, where @f$E@f$ and 
 * @f$D@f$ are the row and column scaling diagonal matrices respectively. The upper and lower bounds 
 * are also scaled with @f$E@f$, thus @f$\bar{b}_\textrm{min}, \bar{b}_\textrm{max} \leftarrow 
 * E b_\textrm{min}, E b_\textrm{max}@f$. The primal variables are transformed using @f$D^{-1}@f$, 
 * resulting in @f$\bar{x}\leftarrow D^{-1}x@f$. Therefore, also the cost matrix Q and vector q have 
 * to be scaled with @f$D@f$, @f$\bar{Q}\leftarrow DQD@f$ and @f$\bar{q}\leftarrow Dq@f$. Finally the 
 * objective function is scaled with a scalar @f$c@f$, thus @f$\bar{Q}, \bar{q} \leftarrow c\bar{Q}, 
 * c\bar{q} @f$, where @f$c=1/\textrm{max}(1, \|D(Qx^0+q)\|_\infty)@f$. The dual variables in the scaled problem
 * become @f$\bar{y} \leftarrow c E^{-1} y@f$. The diagonals of @f$D@f$ and @f$E@f$ are stored and used in the
 * remainder of the problem.
 * 
 *  
 * 
 * @param   work Workspace
 */
void scale_data(QPALMWorkspace *work);

/**
 * Unscale the problem data. See scale_data.
 * 
 * @param   work Workspace
 */
void unscale_data(QPALMWorkspace *work);

# ifdef __cplusplus
}
# endif 

#endif // ifndef SCALING_H
