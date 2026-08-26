! THIS VERSION: 25/04/2022 AT 13:45 GMT
! Nick Gould (nick.gould@stfc.ac.uk)

!> @defgroup qpalm-fortran-grp Fortran Interface
!! This is the Fortran interface of the QPALM solver.

!>  Fortran interface to the C package QPALM, with the aim to
!>  minimize the objective function
!>
!>      1/2 x' H x + g' x + f
!>
!>  subject to the constraints
!>
!>      cl <= A x <= cu,
!>
!>  where any of the bounds cl, cu may be infinite, See the comments to
!>  the interface block, qpalm_fortran, below for details on how to call the
!>  subroutine, and check qpalm_fortran_example.f90 for an example of use
!>
!>  @ingroup qpalm-fortran-grp

MODULE QPALM_fiface

   USE iso_c_binding, ONLY : C_FLOAT, C_DOUBLE, C_CHAR, C_INT32_T, C_INT64_T

   IMPLICIT NONE

   PUBLIC

   !-------------------------------------------
   !   I n r e g e r  a n d  r e a l  k i n d s
   !-------------------------------------------

   !  integer and real kinds for problem data

#ifdef QPALM_FORTRAN_64BIT_INDICES
   INTEGER, PARAMETER :: integer_kind = C_INT64_T
#else
   INTEGER, PARAMETER :: integer_kind = C_INT32_T
#endif

#ifdef QPALM_FORTRAN_SINGLE_PRECISION
   INTEGER, PARAMETER :: real_kind = C_FLOAT
#else
   INTEGER, PARAMETER :: real_kind = C_DOUBLE
#endif

   !  integer and real kinds for qpalm-related data

#ifdef LADEL_64BIT_INDICES
   INTEGER, PARAMETER :: integer_kind_qpalm = C_INT64_T
#else
   INTEGER, PARAMETER :: integer_kind_qpalm = C_INT32_T
#endif

#ifdef LADEL_SINGLE_PRECISION
   INTEGER, PARAMETER :: real_kind_qpalm = C_FLOAT
#else
   INTEGER, PARAMETER :: real_kind_qpalm = C_DOUBLE
#endif

   !----------------------
   !   P a r a m e t e r s
   !----------------------

   REAL ( KIND = real_kind_qpalm ), PARAMETER :: ten = 10.0_real_kind_qpalm

   !-------------------------------------------------
   !  D e r i v e d   t y p e   d e f i n i t i o n s
   !-------------------------------------------------

   !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !> settings derived type with component defaults bound to C’s QPALMSettings
   !> @see @ref QPALMSettings
   !> @ingroup qpalm-fortran-grp
   !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   TYPE, BIND( C ), PUBLIC :: QPALM_settings

      !>  maximum number of iterations: > 0

      INTEGER ( KIND = integer_kind_qpalm ) :: max_iter = 10000

      !>  maximum number of iterations per subproblem: > 0

      INTEGER ( KIND = integer_kind_qpalm ) :: inner_max_iter = 100

      !>  absolute convergence tolerance: >= 0,
      !>   either eps_abs or eps_rel must be > 0

      REAL ( KIND = real_kind_qpalm ) :: eps_abs = ten ** ( - 4 )

      !>  relative convergence tolerance: >= 0,
      !>   either eps_abs or eps_rel must be > 0

      REAL ( KIND = real_kind_qpalm ) :: eps_rel = ten ** ( - 4 )

      !>  intermediate absolute convergence tolerance: >= 0,
      !>   either eps_abs_in or eps_rel_in must be > 0

      REAL ( KIND = real_kind_qpalm ) :: eps_abs_in = 1.0_real_kind_qpalm

      !>  intermediate relative convergence tolerance: >= 0,
      !>   either eps_abs_in or eps_rel_in must be > 0

      REAL ( KIND = real_kind_qpalm ) :: eps_rel_in = 1.0_real_kind_qpalm

      !>  tolerance scaling factor: 0 < rho < 1

      REAL ( KIND = real_kind_qpalm ) :: rho = 0.1_real_kind_qpalm

      !>  primal infeasibility tolerance: >= 0

      REAL ( KIND = real_kind_qpalm ) :: eps_prim_inf = ten ** ( - 5 )

      !>  dual infeasibility tolerance: >= 0

      REAL ( KIND = real_kind_qpalm ) :: eps_dual_inf = ten ** ( - 5 )

      !>  penalty update criterion parameter: <= 1

      REAL ( KIND = real_kind_qpalm ) :: theta = 0.25_real_kind_qpalm

      !>  penalty update factor: > 1

      REAL ( KIND = real_kind_qpalm ) :: delta = 100.0_real_kind_qpalm

      !>  penalty factor cap: > 0

      REAL ( KIND = real_kind_qpalm ) :: sigma_max = ten ** 9

      !>  initial penalty parameter (guideline): > 0

      REAL ( KIND = real_kind_qpalm ) :: sigma_init = 20.0_real_kind_qpalm

      !>  boolean, use proximal method of multipliers or not: in {0,1}

      INTEGER ( KIND = integer_kind_qpalm ) :: proximal = 1

      !>  initial proximal penalty parameter: > 0

      REAL ( KIND = real_kind_qpalm ) :: gamma_init = ten ** 7

      !>  proximal penalty update factor: >= 1

      REAL ( KIND = real_kind_qpalm ) :: gamma_upd = 10.0_real_kind_qpalm

      !>  proximal penalty parameter cap: >= gamma_init

      REAL ( KIND = real_kind_qpalm ) :: gamma_max = ten ** 7

      !>  scaling iterations, if 0 then scaling is disabled: >= 0

      INTEGER ( KIND = integer_kind_qpalm ) :: scaling = 10

      !>  boolean, indicates whether the QP is nonconvex: in {0,1}

      INTEGER ( KIND = integer_kind_qpalm ) :: nonconvex = 0

      !>  boolean, write out progress: in {0,1}

      INTEGER ( KIND = integer_kind_qpalm ) :: verbose = 1

      !>  frequency of printing: > 0

      INTEGER ( KIND = integer_kind_qpalm ) :: print_iter = 1

      !>  boolean, warm start: in {0,1}

      INTEGER ( KIND = integer_kind_qpalm ) :: warm_start = 0

      !>  frequency of performing a complete Cholesky factorization: > 0

      INTEGER ( KIND = integer_kind_qpalm ) :: reset_newton_iter = 10000

      !>  boolean, enable termination based on dual objective (useful
      !>   in branch and bound): in {0,1}

      INTEGER ( KIND = integer_kind_qpalm ) :: enable_dual_termination = 0

      !>  termination value for the dual objective (useful in branch and bound)

      REAL ( KIND = real_kind_qpalm ) :: dual_objective_limit = ten ** 20

      !>  time limit: > 0

      REAL ( KIND = real_kind_qpalm ) :: time_limit = ten ** 20

      !>  ordering method for factorization:
      !>   0  No ordering is performed during the symbolic part of the factorization
      !>   1  Ordering method during the symbolic part of the factorization
      !>   2  The ordering was computed previously and is already stored

      INTEGER ( KIND = integer_kind_qpalm ) :: ordering = 1

      !>  factorize KKT or Schur complement:
      !>    0  factorize the kkt system
      !>    1  factorize the Schur complement
      !>    2  select automatically between kkt system and schur complemen

      INTEGER ( KIND = integer_kind_qpalm ) :: factorization_method = 2

      !>  maximum rank for the sparse factorization update

      INTEGER ( KIND = integer_kind_qpalm ) :: max_rank_update = 160

      !>  maximum rank (relative to n+m) for the factorization update

      REAL ( KIND = real_kind_qpalm ) ::                                     &
         max_rank_update_fraction = 0.1_real_kind_qpalm

   END TYPE QPALM_settings

   !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !>   info derived type with component defaults bound to C’s QPALMInfo
   !>   @see @ref QPALMInfo
   !>   @ingroup qpalm-fortran-grp
   !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   TYPE, BIND( C ), PUBLIC :: QPALM_info

      !>  number of iterations taken

      INTEGER ( KIND = integer_kind_qpalm ) :: iter

      !>  number of outer iterations (i.e. dual updates)

      INTEGER ( KIND = integer_kind_qpalm ) :: iter_out

      !>  status string, e.g. 'solved'

      CHARACTER ( KIND = C_CHAR ), DIMENSION( 31 ) :: status

      !>  status as integer:
      !>  *   1 the problem is solved to optimality given the specified tolerances
      !>  *   2 the problem has a dual objective that is higher than the specified bound
      !>  *   0 an error has occured (this error should automatically be printed)
      !>  *  -2 termination due to reaching the maximum number of iterations
      !>  *  -3 the problem is primal infeasible
      !>  *  -4 the problem is dual infeasible
      !>  *  -5 the problem’s runtime has exceeded the specified time limit
      !>  * -10 the problem is unsolved. Only setup function has been called

      INTEGER ( KIND = integer_kind_qpalm ) :: status_val

      !>  norm of primal residual

      REAL ( KIND = real_kind_qpalm ) :: pri_res_norm

      !>  norm of dual residual

      REAL ( KIND = real_kind_qpalm ) :: dua_res_norm

      !>  norm of intermediate dual residual (minus proximal term)

      REAL ( KIND = real_kind_qpalm ) :: dua2_res_norm

      !>  objective function value

      REAL ( KIND = real_kind_qpalm ) :: objective

      !>  dual objective function value (= NaN if enable_dual_termination is false)

      REAL ( KIND = real_kind_qpalm ) :: dual_objective

      !>  time taken for setup phase (seconds)

      REAL ( KIND = real_kind_qpalm ) :: setup_time

      !>  time taken for solve phase (seconds)

      REAL ( KIND = real_kind_qpalm ) :: solve_time

      !>  total time (seconds)

      REAL ( KIND = real_kind_qpalm ) :: run_time

   END TYPE QPALM_info

   !---------------------------------
   !   I n t e r f a c e  B l o c k s
   !---------------------------------

   INTERFACE
      !> Invoke the QPALM solver.
      !> @see @ref qpalm_solve
      !> @ingroup qpalm-fortran-grp
      SUBROUTINE qpalm_fortran( n, m, hne, hptr, hrow, hval, g, f,             &
         ane, aptr, arow, aval, cl, cu, settings,       &
         x, y, info ) BIND( C, NAME = 'qpalm_fortran_c' )

         !    dummy arguments

         IMPORT :: integer_kind, real_kind, QPALM_settings, QPALM_info

         !>  the number of variables

         INTEGER ( KIND = integer_kind ), INTENT( IN ), VALUE :: n

         !>  the number of constraints

         INTEGER ( KIND = integer_kind ), INTENT( IN ), VALUE :: m

         !>  the number of nonzeros in the upper triangular part of the objective
         !>  Hessian, H

         INTEGER ( KIND = integer_kind ), INTENT( IN ), VALUE ::  hne

         !>  the (1-based) row indices of the upper triangular part of H when H is
         !>  stored by columns. Columns are stored consecutively, with column i directly
         !>  before column i+1, for i = 1,...,n-1

         INTEGER ( KIND = integer_kind ), INTENT( IN ), DIMENSION( hne ) :: hrow

         !>  the (1-based) pointers to the first entry in each column of the upper
         !>  triangular part of H, i = 1,...,n, as well as the pointer to one position
         !>  beyond the last entry, stored in  hptr(n+1)

         INTEGER ( KIND = integer_kind ), INTENT( IN ), DIMENSION( n + 1 ) :: hptr

         !>  the values of the nonzeros in the upper triangular part of H, in the same
         !>  order as the row indices stored in hrow

         REAL ( KIND = real_kind ), INTENT( IN ), DIMENSION( hne ) :: hval

         !>  the vector of values of the linear term, g, in the objective

         REAL ( KIND = real_kind ), INTENT( IN ), DIMENSION( n ) :: g

         !>  the value of the constant term, f, in the objective

         REAL ( KIND = real_kind ), INTENT( IN ), VALUE :: f

         !>  the number of nonzeros in the constraint Jacobian, A

         INTEGER ( KIND = integer_kind ), INTENT( IN ), VALUE ::  ane

         !>  the (1-based) row indices of the A when A is stored by columns. Columns
         !>  are stored consecutively, with column i directly before column i+1, for
         !>  i = 1,...,n-1

         INTEGER ( KIND = integer_kind ), INTENT( IN ), DIMENSION( ane ) :: arow

         !>  the (1-based) pointers to the first entry in each column of A, i = 1,...,n,
         !>  as well as the pointer to one position beyond the last entry, stored in
         !>  aptr(n+1)

         INTEGER ( KIND = integer_kind ), INTENT( IN ), DIMENSION( n + 1 ) :: aptr

         !>  the values of the nonzeros in A, in the same order as the row indices
         !>  stored in arow

         REAL ( KIND = real_kind ), INTENT( IN ), DIMENSION( ane ) :: aval

         !>  the vector of lower constraint bounds, cl. An infinite bound should
         !>  be given a value no larger than - QPALM_INFTY = -10^20

         REAL ( KIND = real_kind ), INTENT( IN ), DIMENSION( m ) :: cl

         !>  the vector of upper constraint bounds, cu. An infinite bound should
         !>  be given a value no smaller than QPALM_INFTY = 10^20

         REAL ( KIND = real_kind ), INTENT( IN ), DIMENSION( m ) :: cu

         !>  parameters that are used to control the optimization. See QPALM_settings
         !>  above for details

         TYPE ( QPALM_settings ), INTENT( IN ), VALUE :: settings

         !>  the values of the best primal variables, x, on successful termination

         REAL ( KIND = real_kind ), INTENT( OUT ), DIMENSION( n ) :: x

         !>  the values of the best dual variables, y, on successful termination

         REAL ( KIND = real_kind ), INTENT( OUT ), DIMENSION( m ) :: y

         !>  output information after the optimization. See QPALM_info above for details

         TYPE ( QPALM_info ), INTENT( OUT ) :: info

      END SUBROUTINE qpalm_fortran
   END INTERFACE

END MODULE QPALM_fiface
