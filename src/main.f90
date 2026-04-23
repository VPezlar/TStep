! =============================================================================
! main  --  program entry point for the TStep time-stepping Arnoldi stability
!           analysis (Mathias & Medeiros, 2022).
!
! Pipeline (one invocation drives the full workflow):
!
!   Initialization
!     1. configurationRead        -- parse inputs/inputs.in.
!     2. validate_paths           -- check baseflow_grid, baseflow_field,
!                                     stability_dir.
!     3. read_flowfield           -- load base flow from the SACRED archive:
!                                     grid from baseflow_grid, fields (rho,
!                                     p, T, U, V, W) from baseflow_field.
!     4. seed_stability_initial   -- mkdir <stab>/1/ and copy the four
!                                     field files from baseflow_field so
!                                     write_flowfield has header templates.
!     5. run_simulation (warmup)  -- advance <stab>/1/ -> <stab>/<1+TTime>/
!                                     to digest the raw snapshot into a
!                                     discretely-consistent flow field.
!     6. promote_to_initial_state -- cp <stab>/<1+TTime>/{p,rho,T,U} back
!                                     over <stab>/1/. Now /1/ holds q0,
!                                     the Arnoldi base state.
!     7. read_flowfield(/1/)      -- load q0 into memory arrays rho0..W0
!                                     (alive for the whole Arnoldi loop).
!     8. write_flowfield_data     -- dump q0 to <stability_dir>/flowfield.csv
!                                     for external inspection.
!     9. IF frechet_order == 1:
!            run_simulation       -- second CFD call produces F(q0) at
!                                     <stab>/<1+TTime>/.
!            read_flowfield(/1+T/)-- pack F(q0) into F_q0_vec so the matvec
!                                     can skip the alpha=0 CFD call on
!                                     every Arnoldi iteration.
!        ELSE (even orders): F_q0_vec stays zero; matvec never touches it.
!    10. initial_disturbance      -- random unit-norm v_1 scaled by eps_0.
!
!   Arnoldi loop (inside arnoldi_eigenvalues; krylov_size iterations)
!     For each j = 1..krylov_size:
!       apply_linearized_operator(v_j, w)  -- sweeps Frechet stencil nodes;
!                                             each non-zero node triggers a
!                                             write -> CFD run -> read.
!       modified Gram-Schmidt + reorth      -- builds H(:,j), v_{j+1}.
!     Finalize: DGEEV on H_m, lift Ritz vectors = V * y, sort.
!
!   Output
!    11. write_eigen_files         -- eigenvalues.dat and eigenvectors.dat
!                                      in <stability_dir>/.
!    12. cleanup_allocations       -- guarded DEALLOCATE + STOP 0.
!
! SACRED: baseflow_grid and baseflow_field are read-only. TStep never
!          writes into them. Everything else lives under stability_dir.
! =============================================================================
PROGRAM main
    USE accuracy
    USE setup
    USE read_flow
    USE write_flow
    USE write_output
    USE write_eigendata
    USE call_CFD
    USE random_disturbance
    USE error_handling
    USE variables
    USE state_vector
    USE Arnoldi
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    ! C interface for setenv (used by set_blas_threads utility below)
    INTERFACE
        FUNCTION c_setenv(name, value, overwrite) BIND(C, NAME="setenv")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: c_setenv
            CHARACTER(KIND=C_CHAR), DIMENSION(*) :: name, value
            INTEGER(C_INT), VALUE :: overwrite
        END FUNCTION c_setenv
    END INTERFACE

    ! --- Raw baseflow (read from the SACRED archive; kept for grid + CSV) ---
    REAL(rk), DIMENSION(:), ALLOCATABLE :: rho_raw, p_raw, T_raw, U_raw, V_raw, W_raw
    REAL(rk), DIMENSION(:), ALLOCATABLE :: Xgrid, Ygrid, Zgrid

    ! --- Base state q0 after the CFD warmup pass (Arnoldi's Frechet origin)
    REAL(rk), DIMENSION(:), ALLOCATABLE :: rho0, p0, T0, U0, V0, W0

    ! --- F(q0) cached in memory for the frechet_order=1 alpha=0 shortcut ---
    REAL(rk), DIMENSION(:), ALLOCATABLE :: F_q0_rho, F_q0_p, F_q0_T
    REAL(rk), DIMENSION(:), ALLOCATABLE :: F_q0_U, F_q0_V, F_q0_W
    REAL(rk), DIMENSION(:), ALLOCATABLE :: F_q0_vec
    LOGICAL                             :: has_F_q0

    ! --- Initial Krylov vector v_1 (random disturbance, length NVARS*Nmesh) --
    REAL(rk), DIMENSION(:), ALLOCATABLE :: v1

    ! --- Arnoldi outputs ---
    COMPLEX(rk), DIMENSION(:),   ALLOCATABLE :: eigenvalues
    COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: eigenvectors

    INTEGER(ik) :: Nmesh, dc_tmp, vlen
    INTEGER(ik) :: error_status

    ! =========================================================================
    ! INITIALIZATION
    ! =========================================================================

    ! --- 1. Read configuration ---
    CALL configurationRead(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_CONFIG)
        STOP ERR_MAIN_CONFIG
    END IF

    ! --- 2. Validate raw path inputs ---
    CALL validate_paths(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_CONFIG)
        STOP ERR_MAIN_CONFIG
    END IF

    ! --- 3. Read raw baseflow from the SACRED archive ---x
    CALL read_flowfield(rho_raw, p_raw, T_raw, U_raw, V_raw, W_raw, &
                        Xgrid, Ygrid, Zgrid, Nmesh, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_READ_FLOW)
        CALL cleanup_allocations()
        STOP ERR_MAIN_READ_FLOW
    END IF
    WRITE(*,'(A,I0)') '[init] Nmesh = ', Nmesh

    ! --- 4. Seed <stability_dir>/1/ with the four baseflow field files ---
    CALL seed_stability_initial(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! --- 5. CFD warmup pass: <stab>/1/ -> <stab>/<1+TTime>/ ---
    WRITE(*,'(A)') '[init] CFD warmup pass (digest raw baseflow)'
    CALL clear_endpoint_folder(error_status)
    CALL run_simulation(TRIM(COMMAND_RUN), error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_EXT_CMD)
        CALL cleanup_allocations()
        STOP ERR_MAIN_EXT_CMD
    END IF

    ! Verify the solver actually produced an endpoint folder, and pin down
    ! its name via find_endpoint_folder(). This decouples TStep from
    ! OpenFOAM's timePrecision / round-off behaviour: whatever the solver
    ! called its output folder (1.1/, 1.099999999999989/, ...), we pick it
    ! up by scanning <stability_dir> for the latest numeric-named folder.
    BLOCK
        CHARACTER(len=256) :: ep_warmup
        LOGICAL :: ep_file_exists
        CALL find_endpoint_folder(ep_warmup, error_status)
        IF (error_status /= 0) THEN
            CALL cleanup_allocations()
            STOP ERR_SETUP_INVALID_PARAM
        END IF
        INQUIRE(FILE=TRIM(ep_warmup)//'p', EXIST=ep_file_exists)
        IF (.NOT. ep_file_exists) THEN
            CALL log_error(ERR_SETUP_INVALID_PARAM, &
                'Warmup endpoint folder '//TRIM(ep_warmup)//&
                ' exists but has no p file. Check controlDict.')
            CALL cleanup_allocations()
            STOP ERR_SETUP_INVALID_PARAM
        END IF
    END BLOCK

    ! --- 6. Promote warmup output to /1/: now /1/ holds q0 ---
    CALL promote_to_initial_state(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! --- 7. Load q0 into memory (base state for the Arnoldi loop) ---
    CALL read_flowfield(rho0, p0, T0, U0, V0, W0, &
                        Xgrid, Ygrid, Zgrid, dc_tmp, error_status, &
                        path_override=stability_time_dir(TSTEP_INITIAL_TIME))
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_READ_FLOW)
        CALL cleanup_allocations()
        STOP ERR_MAIN_READ_FLOW
    END IF
    IF (dc_tmp /= Nmesh) THEN
        CALL log_error(ERR_MAIN_READ_FLOW, 'q0 read returned different Nmesh than baseflow')
        CALL cleanup_allocations()
        STOP ERR_MAIN_READ_FLOW
    END IF

    ! --- 8. Dump q0 to <stability_dir>/flowfield.csv for inspection ---
    CALL write_flowfield_data(Xgrid, Ygrid, Zgrid, rho0, p0, T0, &
                              U0, V0, W0, Nmesh, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! --- 9. For frechet_order == 1: precompute F(q0) and cache it ---
    vlen = state_length(Nmesh)
    ALLOCATE(F_q0_vec(vlen), STAT=error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT, 'F_q0_vec allocation failed')
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF
    F_q0_vec = 0.0_rk
    has_F_q0 = .FALSE.

    IF (frechet_order == 1) THEN
        WRITE(*,'(A)') '[init] frechet_order=1: computing F(q0) cache'
        CALL clear_endpoint_folder(error_status)
        CALL run_simulation(TRIM(COMMAND_RUN), error_status)
        IF (error_status /= 0) THEN
            CALL log_error(ERR_MAIN_EXT_CMD)
            CALL cleanup_allocations()
            STOP ERR_MAIN_EXT_CMD
        END IF

        BLOCK
            CHARACTER(len=256) :: ep_cache
            LOGICAL :: fq0_exists
            CALL find_endpoint_folder(ep_cache, error_status)
            IF (error_status /= 0) THEN
                CALL cleanup_allocations()
                STOP ERR_SETUP_INVALID_PARAM
            END IF
            INQUIRE(FILE=TRIM(ep_cache)//'p', EXIST=fq0_exists)
            IF (.NOT. fq0_exists) THEN
                CALL log_error(ERR_SETUP_INVALID_PARAM, &
                    'F(q0) cache endpoint folder '//TRIM(ep_cache)//&
                    ' has no p file. Check controlDict.')
                CALL cleanup_allocations()
                STOP ERR_SETUP_INVALID_PARAM
            END IF

            CALL read_flowfield(F_q0_rho, F_q0_p, F_q0_T, F_q0_U, F_q0_V, F_q0_W, &
                                Xgrid, Ygrid, Zgrid, dc_tmp, error_status, &
                                path_override=ep_cache)
        END BLOCK
        IF (error_status /= 0) THEN
            CALL log_error(ERR_MAIN_READ_FLOW)
            CALL cleanup_allocations()
            STOP ERR_MAIN_READ_FLOW
        END IF

        BLOCK
            REAL(rk) :: eps_s, norm_q0
            REAL(rk), ALLOCATABLE :: q0_vec(:)
            ALLOCATE(q0_vec(vlen))
            CALL pack_state(rho0, p0, T0, U0, V0, W0, q0_vec, error_status)
            norm_q0 = NORM2(q0_vec)
            eps_s = NORM2(F_q0_vec - q0_vec)
            WRITE(*,*) 'TEST FOR NOISE STARTS HERE', norm_q0
            WRITE(*,'(A,ES12.4)') '[diag] ||q0||        = ', norm_q0
            WRITE(*,'(A,ES12.4)') '[diag] ||F(q0)-q0|| = ', eps_s
            WRITE(*,'(A,ES12.4)') '[diag] relative εₛ  = ', eps_s/norm_q0
            DEALLOCATE(q0_vec)

        END BLOCK

        CALL pack_state(F_q0_rho, F_q0_p, F_q0_T, F_q0_U, F_q0_V, F_q0_W, &
                        F_q0_vec, error_status)
        IF (error_status /= 0) THEN
            CALL log_error(ERR_MAIN_DISTURBANCE, 'pack_state on F(q0) failed')
            CALL cleanup_allocations()
            STOP ERR_MAIN_DISTURBANCE
        END IF
        has_F_q0 = .TRUE.
    END IF

    ! --- 10. Initial Krylov vector v_1 (random unit-norm, scaled by eps_0) --
    CALL initial_disturbance(vlen, eps_0, v1, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_DISTURBANCE)
        CALL cleanup_allocations()
        STOP ERR_MAIN_DISTURBANCE
    END IF

    ! =========================================================================
    ! ARNOLDI LOOP  (matvec + Gram-Schmidt + DGEEV + sort)
    ! =========================================================================
    WRITE(*,'(A,I0,A)') '[arnoldi] starting Krylov loop (m = ', krylov_size, ')'
    CALL arnoldi_eigenvalues(v1, krylov_size, frechet_order, eps_0, TTime, &
                             rho0, p0, T0, U0, V0, W0, &
                             F_q0_vec, has_F_q0, &
                             eigenvalues, eigenvectors, error_status, &
                             sort_by=eigenvalue_sort_by)
    IF (error_status /= 0) THEN
        CALL cleanup_allocations()
        STOP error_status
    END IF

    ! =========================================================================
    ! OUTPUT
    ! =========================================================================
    CALL write_eigen_files(eigenvalues, eigenvectors, Xgrid, Ygrid, Zgrid, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    CALL cleanup_allocations()

CONTAINS

    SUBROUTINE set_blas_threads(nthreads)
        ! Sets thread count for BLAS/LAPACK via environment variables.
        ! Unused today; intended for tuning the real Arnoldi run once profile
        ! data is available.
        INTEGER(ik), INTENT(IN) :: nthreads
        CHARACTER(len=20)  :: threads_str
        CHARACTER(len=100) :: env_value
        INTEGER(C_INT)     :: result
        INTEGER            :: env_length, env_status

        IF (nthreads <= 0) THEN
            WRITE(*,*)
            WRITE(*,*) 'BLAS threading: AUTO mode (using all available cores)'
        ELSE
            WRITE(threads_str, '(I0)') nthreads
            result = c_setenv('OMP_NUM_THREADS'//C_NULL_CHAR, &
                              TRIM(threads_str)//C_NULL_CHAR, 1_C_INT)
            result = c_setenv('OPENBLAS_NUM_THREADS'//C_NULL_CHAR, &
                              TRIM(threads_str)//C_NULL_CHAR, 1_C_INT)
            result = c_setenv('MKL_NUM_THREADS'//C_NULL_CHAR, &
                              TRIM(threads_str)//C_NULL_CHAR, 1_C_INT)

            WRITE(*,*)
            WRITE(*,'(A,I0,A)') 'BLAS threading: Set to ', nthreads, ' thread(s)'
            CALL GET_ENVIRONMENT_VARIABLE('OMP_NUM_THREADS', env_value, env_length, env_status)
            IF (env_status == 0) THEN
                WRITE(*,'(A,A)') 'Verified OMP_NUM_THREADS = ', TRIM(env_value)
            END IF
        END IF
    END SUBROUTINE set_blas_threads


    SUBROUTINE cleanup_allocations()
        ! Raw baseflow
        IF (ALLOCATED(rho_raw)) DEALLOCATE(rho_raw)
        IF (ALLOCATED(p_raw))   DEALLOCATE(p_raw)
        IF (ALLOCATED(T_raw))   DEALLOCATE(T_raw)
        IF (ALLOCATED(U_raw))   DEALLOCATE(U_raw)
        IF (ALLOCATED(V_raw))   DEALLOCATE(V_raw)
        IF (ALLOCATED(W_raw))   DEALLOCATE(W_raw)
        IF (ALLOCATED(Xgrid))   DEALLOCATE(Xgrid)
        IF (ALLOCATED(Ygrid))   DEALLOCATE(Ygrid)
        IF (ALLOCATED(Zgrid))   DEALLOCATE(Zgrid)

        ! Base state q0
        IF (ALLOCATED(rho0))    DEALLOCATE(rho0)
        IF (ALLOCATED(p0))      DEALLOCATE(p0)
        IF (ALLOCATED(T0))      DEALLOCATE(T0)
        IF (ALLOCATED(U0))      DEALLOCATE(U0)
        IF (ALLOCATED(V0))      DEALLOCATE(V0)
        IF (ALLOCATED(W0))      DEALLOCATE(W0)

        ! F(q0) cache
        IF (ALLOCATED(F_q0_rho)) DEALLOCATE(F_q0_rho)
        IF (ALLOCATED(F_q0_p))   DEALLOCATE(F_q0_p)
        IF (ALLOCATED(F_q0_T))   DEALLOCATE(F_q0_T)
        IF (ALLOCATED(F_q0_U))   DEALLOCATE(F_q0_U)
        IF (ALLOCATED(F_q0_V))   DEALLOCATE(F_q0_V)
        IF (ALLOCATED(F_q0_W))   DEALLOCATE(F_q0_W)
        IF (ALLOCATED(F_q0_vec)) DEALLOCATE(F_q0_vec)

        ! Krylov + Arnoldi outputs
        IF (ALLOCATED(v1))           DEALLOCATE(v1)
        IF (ALLOCATED(eigenvalues))  DEALLOCATE(eigenvalues)
        IF (ALLOCATED(eigenvectors)) DEALLOCATE(eigenvectors)
    END SUBROUTINE cleanup_allocations

END PROGRAM main
