PROGRAM main
    USE accuracy
    USE setup
    USE read_flow
    USE write_flow
    USE write_output
    USE call_CFD
    USE random_disturbance
    USE error_handling
    USE variables
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

    REAL(rk), DIMENSION(:), ALLOCATABLE :: rho_in, p_in, T_in, U_in, V_in, W_in
    REAL(rk), DIMENSION(:), ALLOCATABLE :: Xgrid, Ygrid, Zgrid, pert_0
    INTEGER(ik) :: data_count
    INTEGER(ik) :: error_status, STATUS_CODE
    INTEGER(ik) :: unit_num, i

    ! --- Read configuration ---
    CALL configurationRead(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_CONFIG)
        STOP ERR_MAIN_CONFIG
    END IF

    ! --- Execute external command (run CFD solver once to generate base state) ---
    CALL run_simulation(COMMAND_RUN, STATUS_CODE)
    IF (STATUS_CODE /= 0) THEN
        CALL log_error(ERR_MAIN_EXT_CMD)
        STOP ERR_MAIN_EXT_CMD
    END IF

    WRITE(*,*) 'PROCEEDING TO NEXT STEP.'

    ! --- Read Flowfield ---
    CALL read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, &
                        Xgrid, Ygrid, Zgrid, data_count, error_status)

    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_READ_FLOW)
        CALL cleanup_allocations()
        STOP ERR_MAIN_READ_FLOW
    END IF

    ! --- Generate initial disturbance ---
    CALL initial_disturbance(data_count, dist_mag, pert_0, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_DISTURBANCE)
        CALL cleanup_allocations()
        STOP ERR_MAIN_DISTURBANCE
    END IF

    ! --- TEMPORARY: Write disturbance vector to file for testing ---
    ! TODO: Remove this section after validation is complete
    CALL get_unit(unit_num)
    OPEN(UNIT=unit_num, FILE='disturbance.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'WARNING: Failed to open disturbance.txt for writing (IOSTAT:', error_status, ')'
        WRITE(*,*) 'Continuing without writing disturbance file...'
    ELSE
        ! Write vector size first (helpful for later reading)
        WRITE(unit_num, *) SIZE(pert_0)
        ! Write all elements of the vector, one per line
        DO i = 1, data_count
            WRITE(unit_num, *) pert_0(i)
        END DO
        CLOSE(unit_num)
        WRITE(*,*) 'SUCCESS: Wrote disturbance vector to disturbance.txt.'
    END IF
    ! --- End TEMPORARY section ---

    ! --- Write perturbed flowfield (base + disturbance) ---
    CALL write_flowfield(rho_in + pert_0, &
                         p_in   + pert_0, &
                         T_in   + pert_0, &
                         U_in   + pert_0, &
                         V_in   + pert_0, &
                         W_in   + pert_0, &
                         data_count, error_status)

    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! --- Write flowfield data to CSV ---
    CALL write_flowfield_data(Xgrid, Ygrid, Zgrid, rho_in, p_in, T_in, &
                              U_in, V_in, W_in, data_count, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! -------------------------------------------------------------------------
    ! NOTE: The Arnoldi eigenvalue validation test (n=7500 reference problem)
    ! has been removed. The real Arnoldi call needs a working Frechet-derivative
    ! matvec in Arnoldi.f90 (currently a stub). Re-add the invocation here once
    ! that stub is implemented. See:
    !   - src/Arnoldi.f90 :: apply_linearized_operator (TODO stub)
    !   - src/call_CFD.f90 :: run_simulation (solver advance)
    !   - ERR_ARNOLDI_NOT_IMPLEMENTED in error_handling.f90
    ! -------------------------------------------------------------------------

    ! --- Cleanup allocated memory ---
    CALL cleanup_allocations()

CONTAINS

    SUBROUTINE set_blas_threads(nthreads)
        ! Sets thread count for BLAS/LAPACK via environment variables.
        ! Unused today; intended for the real Arnoldi run once the matvec lands.
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
        IF (ALLOCATED(p_in))   DEALLOCATE(p_in)
        IF (ALLOCATED(rho_in)) DEALLOCATE(rho_in)
        IF (ALLOCATED(T_in))   DEALLOCATE(T_in)
        IF (ALLOCATED(U_in))   DEALLOCATE(U_in)
        IF (ALLOCATED(V_in))   DEALLOCATE(V_in)
        IF (ALLOCATED(W_in))   DEALLOCATE(W_in)
        IF (ALLOCATED(Xgrid))  DEALLOCATE(Xgrid)
        IF (ALLOCATED(Ygrid))  DEALLOCATE(Ygrid)
        IF (ALLOCATED(Zgrid))  DEALLOCATE(Zgrid)
        IF (ALLOCATED(pert_0)) DEALLOCATE(pert_0)
    END SUBROUTINE cleanup_allocations

END PROGRAM main
