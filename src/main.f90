PROGRAM main
    USE accuracy
    USE setup
    USE read_flow
    USE write_flow
    USE write_output
    USE call_CFD
    USE random_disturbance
    USE error_handling
    USE Arnoldi
    USE write_eigendata
    USE variables

    IMPLICIT NONE

    REAL(rk), DIMENSION(:), ALLOCATABLE :: rho_in, p_in, T_in, U_in, V_in, W_in, Xgrid, Ygrid, Zgrid, pert_0
    INTEGER(ik) :: data_count
    INTEGER(ik) :: error_status, STATUS_CODE
    INTEGER(ik) :: unit_num, i
    
    ! Arnoldi variables
    COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: eigenvalues
    COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: eigenvectors
    REAL(rk), DIMENSION(:), ALLOCATABLE :: v_normalized

    ! Read configuration
    CALL configurationRead(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_CONFIG)
        STOP ERR_MAIN_CONFIG
    END IF

    ! Execute external command
    CALL run_simulation(COMMAND_RUN, STATUS_CODE)
    IF (STATUS_CODE /= 0) THEN
        CALL log_error(ERR_MAIN_EXT_CMD)
        STOP ERR_MAIN_EXT_CMD
    END IF

    ! Continue with other code if STATUS_CODE is 0
    WRITE(*,*) 'PROCEEDING TO NEXT STEP.'

    ! --- Read Flowfield ---
    ! Read flowfield data
    CALL read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, &
                        Xgrid, Ygrid, Zgrid, data_count, error_status)
    
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_READ_FLOW)
        CALL cleanup_allocations()
        STOP ERR_MAIN_READ_FLOW
    END IF

    ! Generate initial disturbance
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

    ! --- Write Flowfield ---
    CALL write_flowfield(rho_in + pert_0, &
                         p_in + pert_0, &
                         T_in + pert_0, &
                         U_in + pert_0, &
                         V_in + pert_0, &
                         W_in + pert_0, &
                         data_count, error_status)

    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! Write flowfield data to file
    CALL write_flowfield_data(Xgrid, Ygrid, Zgrid, rho_in, p_in, T_in, &
                              U_in, V_in, W_in, data_count, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
    END IF

    ! ===================================================================
    ! --- ARNOLDI EIGENVALUE COMPUTATION (REFERENCE TEST) ---
    ! ===================================================================
    ! This test EXACTLY replicates the reference implementation:
    ! https://relate.cs.illinois.edu/.../Arnoldi%20iteration.html
    ! 
    ! Test matrix: n=25, eigenvalues = [1, 2, 3, ..., 25]
    ! Construction: A = eigvecs @ diag(eigvals) @ inv(eigvecs)
    ! Expected: Arnoldi recovers eigenvalues 25, 24, 23, ..., 1
    ! ===================================================================
    
    WRITE(*,*)
    WRITE(*,*) '======================================================='
    WRITE(*,*) 'ARNOLDI TEST - REFERENCE IMPLEMENTATION REPLICATION'
    WRITE(*,*) '======================================================='
    WRITE(*,*)
    WRITE(*,*) 'Reference: CS450 Arnoldi iteration demo'
    WRITE(*,*) 'https://relate.cs.illinois.edu/.../Arnoldi%20iteration.html'
    WRITE(*,*)
    WRITE(*,*) 'Test matrix:'
    WRITE(*,*) '  - Dimension: n = 25'
    WRITE(*,*) '  - Eigenvalues: 1, 2, 3, ..., 25'
    WRITE(*,*) '  - Construction: A = eigvecs @ diag(eigvals) @ inv(eigvecs)'
    WRITE(*,*) '  - Random eigenvector matrix'
    WRITE(*,*)
    WRITE(*,*) 'Expected Ritz values: 25, 24, 23, 22, 21, ..., 1'
    WRITE(*,*) '(Arnoldi finds largest eigenvalues first)'
    WRITE(*,*)
    
    ! HARDCODED test size (matching reference)
    ALLOCATE(v_normalized(25), STAT=error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'ERROR: Failed to allocate v_normalized'
        CALL cleanup_allocations()
        STOP 1
    END IF
    
    ! Initialize random starting vector
    CALL RANDOM_NUMBER(v_normalized)
    v_normalized = v_normalized / NORM2(v_normalized)
    
    WRITE(*,*) '======================================================='
    WRITE(*,*) 'RUNNING ARNOLDI ITERATION'
    WRITE(*,*) '======================================================='
    WRITE(*,*)
    
    ! Call Arnoldi with n=25, m=25
    CALL arnoldi_eigenvalues(v_normalized, 25, frechet_order, eps_0, TTime, &
                            eigenvalues, eigenvectors, error_status, &
                            skip_normalization=.TRUE., sort_by='magnitude')
    
    IF (error_status /= 0) THEN
        CALL log_error(error_status)
        CALL cleanup_allocations()
        STOP error_status
    END IF
    
    ! Display results
    WRITE(*,*)
    WRITE(*,*) '======================================================='
    WRITE(*,*) 'RESULTS'
    WRITE(*,*) '======================================================='
    WRITE(*,*)
    WRITE(*,*) 'Ritz eigenvalues (sorted by magnitude):'
    WRITE(*,*) '-------------------------------------------------------'
    WRITE(*,*) '  #    Real Part      Imag Part      |λ|         Expected'
    WRITE(*,*) '-------------------------------------------------------'
    DO i = 1, 25
        WRITE(*,'(I3,2X,F12.6,2X,F12.6,2X,F12.6,2X,I5)') i, &
            REAL(eigenvalues(i)), AIMAG(eigenvalues(i)), ABS(eigenvalues(i)), 25 - i + 1
    END DO
    WRITE(*,*) '-------------------------------------------------------'
    WRITE(*,*)
    
    ! Write eigenvalues and eigenvectors to files in ../output/
    CALL write_eigen_files(eigenvalues, eigenvectors, error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'WARNING: Failed to write eigendata files'
    END IF
    WRITE(*,*)
    
    ! ===================================================================
    ! --- VALIDATION: Check eigenvalues against analytical values ---
    ! ===================================================================
    CALL validate_eigenvalues()
    
    ! ===================================================================
    ! --- END ARNOLDI SECTION ---
    ! ===================================================================

    ! Cleanup allocated memory
    CALL cleanup_allocations()

CONTAINS

    SUBROUTINE validate_eigenvalues()
        ! Validates against reference test: n=25, eigenvalues 1-25
        REAL(rk) :: max_imag, max_rel_error, rel_error
        REAL(rk) :: expected_eval, computed_val
        INTEGER(ik) :: k
        LOGICAL :: all_real, values_correct
        
        WRITE(*,*) '======================================================='
        WRITE(*,*) 'VALIDATION'
        WRITE(*,*) '======================================================='
        WRITE(*,*)
        
        ! Check 1: Are eigenvalues real?
        max_imag = 0.0_rk
        DO k = 1, 25
            max_imag = MAX(max_imag, ABS(AIMAG(eigenvalues(k))))
        END DO
        
        all_real = (max_imag < 1.0E-10_rk)
        WRITE(*,'(A,ES12.4)') 'Maximum imaginary part: ', max_imag
        IF (all_real) THEN
            WRITE(*,*) 'PASS: All eigenvalues are real (Im(λ) < 1E-10)'
        ELSE
            WRITE(*,*) 'FAIL: Eigenvalues have imaginary parts > 1E-10'
        END IF
        WRITE(*,*)
        
        ! Check 2: Are eigenvalues correct? (25, 24, 23, ..., 1)
        max_rel_error = 0.0_rk
        DO k = 1, 25
            expected_eval = REAL(25 - k + 1, rk)  ! 25, 24, 23, ..., 1
            computed_val = ABS(eigenvalues(k))
            rel_error = ABS(computed_val - expected_eval) / expected_eval * 100.0_rk
            max_rel_error = MAX(max_rel_error, rel_error)
        END DO
        
        values_correct = (max_rel_error < 1.0_rk)
        WRITE(*,'(A,F10.6,A)') 'Maximum relative error: ', max_rel_error, '%'
        IF (max_rel_error < 0.01_rk) THEN
            WRITE(*,*) 'EXCELLENT: Eigenvalues match with < 0.01% error'
        ELSE IF (max_rel_error < 0.1_rk) THEN
            WRITE(*,*) 'VERY GOOD: Eigenvalues match with < 0.1% error'
        ELSE IF (max_rel_error < 1.0_rk) THEN
            WRITE(*,*) 'GOOD: Eigenvalues match with < 1% error'
        ELSE
            WRITE(*,*) 'FAIL: Eigenvalues have > 1% error'
        END IF
        WRITE(*,*)
        
        ! Final verdict
        IF (all_real .AND. values_correct) THEN
            WRITE(*,*) '======================================================='
            WRITE(*,*) 'SUCCESS: ARNOLDI IMPLEMENTATION IS CORRECT!'
            WRITE(*,*) '======================================================='
            WRITE(*,*)
            WRITE(*,*) 'The implementation matches the reference perfectly.'
            WRITE(*,*) 'All eigenvalues recovered correctly.'
        ELSE
            WRITE(*,*) '======================================================='
            WRITE(*,*) 'FAILURE: IMPLEMENTATION DOES NOT MATCH REFERENCE'
            WRITE(*,*) '======================================================='
            IF (.NOT. all_real) THEN
                WRITE(*,*) '  Problem: Eigenvalues have imaginary components'
            END IF
            IF (.NOT. values_correct) THEN
                WRITE(*,*) '  Problem: Eigenvalue errors exceed 1%'
            END IF
        END IF
        WRITE(*,*)
        
    END SUBROUTINE validate_eigenvalues
    
    SUBROUTINE cleanup_allocations()
        IF (ALLOCATED(p_in)) DEALLOCATE(p_in)
        IF (ALLOCATED(rho_in)) DEALLOCATE(rho_in)
        IF (ALLOCATED(T_in)) DEALLOCATE(T_in)
        IF (ALLOCATED(U_in)) DEALLOCATE(U_in)
        IF (ALLOCATED(V_in)) DEALLOCATE(V_in)
        IF (ALLOCATED(W_in)) DEALLOCATE(W_in)
        IF (ALLOCATED(Xgrid)) DEALLOCATE(Xgrid)
        IF (ALLOCATED(Ygrid)) DEALLOCATE(Ygrid)
        IF (ALLOCATED(Zgrid)) DEALLOCATE(Zgrid)
        IF (ALLOCATED(pert_0)) DEALLOCATE(pert_0)
        IF (ALLOCATED(v_normalized)) DEALLOCATE(v_normalized)
        IF (ALLOCATED(eigenvalues)) DEALLOCATE(eigenvalues)
        IF (ALLOCATED(eigenvectors)) DEALLOCATE(eigenvectors)
    END SUBROUTINE cleanup_allocations

END PROGRAM main
