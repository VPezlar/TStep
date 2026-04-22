! =============================================================================
! variables  --  shared configuration globals populated from inputs/inputs.in.
!
! This is the central bag of runtime settings that every other module reads
! after `configurationRead` in setup.f90 has parsed the namelists. Keeping
! them here (rather than passing dozens of arguments) lets the solver-agnostic
! pipeline switch formats or parameters without plumbing changes.
!
! Path model (flat, no derived globals):
!   Three raw inputs from the &OpenFOAM namelist are all TStep needs:
!     baseflow_grid  -- full path to the grid cell-centers file (a FILE;
!                       e.g. '<case>/0/C'). Read once, never written.
!     baseflow_field -- directory holding the canonical base-flow fields
!                       p, rho, T, U (must end with '/'). Read-only archive.
!     stability_dir  -- root of the working/analysis case where stability
!                       analysis runs (must end with '/'). TStep seeds
!                       <stability_dir>/1/ with the base flow before running
!                       the CFD solver; the solver advances it into
!                       <stability_dir>/<time_to_str(1.0+TTime)>/.
!
!   Time-folder paths (<stability_dir>/1/ and <stability_dir>/<1+TTime>/) are
!   composed on demand by callers via setup.f90:stability_time_dir(t). They
!   are NOT stored as module globals -- keeping the model flat avoids the
!   PROTECTED/setter indirection that previously confused the path story.
!
!   OpenFOAM controlDict contract inside <stability_dir>/system/:
!     startTime == TSTEP_INITIAL_TIME        (hardcoded to 1.0)
!     endTime   == TSTEP_INITIAL_TIME+TTime  (varies with TTime in inputs.in)
!     timeFormat    general
!     timePrecision 6
!
! Groups:
!   - Project invariants   (TSTEP_INITIAL_TIME)
!   - Simulation driver    (flow_format, COMMAND_RUN)
!   - Raw paths            (baseflow_grid, baseflow_field, stability_dir,
!                           N_HEADER_grid/var)
!   - Arnoldi parameters   (krylov_size, frechet_order, eigenvalue_sort_by,
!                           eps_0, TTime)
!   - Performance          (num_threads)
! =============================================================================
MODULE variables

    USE accuracy

    IMPLICIT NONE

    ! --- Project-wide invariants --------------------------------------------
    ! TStep always deposits the seeded initial state in <stability_dir>/1/.
    ! The integer '1' here is time_to_str(TSTEP_INITIAL_TIME) under the
    ! OpenFOAM 'timeFormat general, timePrecision 6' convention.
    REAL(rk), PARAMETER :: TSTEP_INITIAL_TIME = 1.0_rk

    ! --- &General namelist: solver-agnostic settings ------------------------
    CHARACTER(len=256) :: flow_format
    CHARACTER(len=256) :: COMMAND_RUN
    REAL(rk)           :: eps_0            ! Frechet perturbation magnitude
                                           ! (replaces the retired dist_mag;
                                           !  also serves as the initial-
                                           !  disturbance L2 norm scale, which
                                           !  is washed out by Arnoldi's own
                                           !  unit-norm renormalization)

    ! --- &OpenFOAM namelist: raw path inputs --------------------------------
    INTEGER(ik)        :: N_HEADER_grid, N_HEADER_var
    CHARACTER(len=256) :: baseflow_grid    ! full path to grid C file (a FILE)
    CHARACTER(len=256) :: baseflow_field   ! dir with p/rho/T/U (ends in '/')
    CHARACTER(len=256) :: stability_dir    ! analysis/working case root ('/')

    ! --- &Arnoldi namelist: algorithm parameters ----------------------------
    INTEGER(ik)       :: krylov_size         ! Krylov subspace size
    INTEGER(ik)       :: frechet_order       ! Order of the frechet derivative
    INTEGER(ik)       :: N_eig_write         ! Number of vectors to be outputted
    CHARACTER(len=20) :: eigenvalue_sort_by  ! switch variable for sorting
    REAL(rk)          :: TTime               ! integration time tau (Mathias eq. 8)

    ! --- &Arnoldi namelist: performance -------------------------------------
    INTEGER(ik)       :: num_threads      ! threads for BLAS/LAPACK (0 = auto)

END MODULE variables
