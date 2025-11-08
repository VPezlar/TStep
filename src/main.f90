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
    ! --- ARNOLDI EIGENVALUE COMPUTATION (TEST/TEMPLATE) ---
    ! ===================================================================
    ! This section demonstrates how to use the Arnoldi module to compute
    ! eigenvalues. Currently uses a hardcoded test matrix.
    ! TODO: Replace matrix-vector product in Arnoldi.f90 with CFD solver calls
    ! ===================================================================
    
    WRITE(*,*) ''
    WRITE(*,*) '==============================================='
    WRITE(*,*) 'ARNOLDI EIGENVALUE COMPUTATION'
    WRITE(*,*) '==============================================='
    WRITE(*,*) ''
    
    ! Normalize the disturbance vector for Arnoldi (requires ||v|| = 1)
    ALLOCATE(v_normalized(SIZE(pert_0)), STAT=error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'ERROR: Failed to allocate v_normalized'
        CALL cleanup_allocations()
        STOP 1
    END IF
    
    v_normalized = pert_0 / NORM2(pert_0)
    WRITE(*,'(A,I0)') 'System dimension n = ', SIZE(v_normalized)
    WRITE(*,'(A,I0)') 'Krylov size m = ', krylov_size
    WRITE(*,'(A,ES15.6)') 'Initial vector normalized: ||v|| = ', NORM2(v_normalized)
    WRITE(*,*) ''
    
    ! Call Arnoldi eigenvalue routine
    ! NOTE: Currently uses hardcoded test matrix. Step 3 in Arnoldi.f90
    !       will be replaced with CFD solver calls for production.
    WRITE(*,*) 'Running Arnoldi iteration...'
    CALL arnoldi_eigenvalues(v_normalized, krylov_size, frechet_order, eps_0, TTime, &
                            eigenvalues, eigenvectors, error_status, &
                            skip_normalization=.TRUE., sort_by=eigenvalue_sort_by)
    
    IF (error_status /= 0) THEN
        CALL log_error(error_status)
        CALL cleanup_allocations()
        STOP error_status
    END IF
    
    ! Display results
    WRITE(*,*) ''
    WRITE(*,*) '==============================================='
    WRITE(*,*) 'SUCCESS: Eigenvalue computation completed!'
    WRITE(*,*) '==============================================='
    WRITE(*,*) ''
    WRITE(*,*) 'Ritz eigenvalues (sorted by descending Im part):'
    WRITE(*,*) '-----------------------------------------------'
    WRITE(*,*) '  #    Real Part         Imag Part         |λ|'
    WRITE(*,*) '-----------------------------------------------'
    DO i = 1, MIN(10, krylov_size)  ! Display first 10 eigenvalues
        WRITE(*,'(I3,2X,ES15.6,2X,ES15.6,2X,ES15.6)') i, &
            REAL(eigenvalues(i)), AIMAG(eigenvalues(i)), ABS(eigenvalues(i))
    END DO
    WRITE(*,*) '-----------------------------------------------'
    WRITE(*,*) ''
    
    ! Write eigenvalues to file
    CALL get_unit(unit_num)
    OPEN(UNIT=unit_num, FILE='eigenvalues.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'WARNING: Failed to open eigenvalues.txt for writing'
    ELSE
        WRITE(unit_num, '(A)') '# Ritz Eigenvalues from Arnoldi Iteration'
        WRITE(unit_num, '(A)') '# Index, Real Part, Imaginary Part, Magnitude'
        DO i = 1, krylov_size
            WRITE(unit_num, '(I5,3ES25.15)') i, REAL(eigenvalues(i)), &
                AIMAG(eigenvalues(i)), ABS(eigenvalues(i))
        END DO
        CLOSE(unit_num)
        WRITE(*,*) 'SUCCESS: Wrote eigenvalues to eigenvalues.txt'
    END IF
    WRITE(*,*) ''
    
    ! ===================================================================
    ! --- VALIDATION: Check eigenvalues against analytical values ---
    ! ===================================================================
    ! The test matrix is tridiagonal with A(i,i)=0, A(i,i±1)=1
    ! Analytical eigenvalues: λ_k = 2*cos(k*π/(n+1)) for k=1,...,n
    ! These are REAL eigenvalues, so Im(λ) should be ~0
    ! ===================================================================
    
    CALL validate_eigenvalues()
    
    ! ===================================================================
    ! --- END ARNOLDI SECTION ---
    ! ===================================================================

    ! Cleanup allocated memory
    CALL cleanup_allocations()

CONTAINS

    SUBROUTINE validate_eigenvalues()
        ! Validates computed Ritz eigenvalues against analytical eigenvalues
        ! of the tridiagonal test matrix
        
        ! --- TEST VALIDATION VARIABLES (delete when removing test matrix) ---
        REAL(rk), DIMENSION(:), ALLOCATABLE :: analytical_evals
        REAL(rk) :: pi
        INTEGER(ik) :: n_analytical
        ! --- END TEST VALIDATION VARIABLES ---
        
        ! Regular variables
        REAL(rk) :: max_error, max_imag, rel_error, max_rel_error
        REAL(rk) :: computed_val, error_abs
        INTEGER(ik) :: k
        LOGICAL :: all_real, validation_passed
        
        WRITE(*,*) '==============================================='
        WRITE(*,*) 'EIGENVALUE VALIDATION'
        WRITE(*,*) '==============================================='
        WRITE(*,*) ''
        
        n_analytical = SIZE(v_normalized)
        
        ! Compute analytical eigenvalues for diagonal test matrix: λ_k = n, n-1, ..., 2, 1
        ALLOCATE(analytical_evals(krylov_size))
        DO k = 1, krylov_size
            analytical_evals(k) = REAL(n_analytical - k + 1, rk)
        END DO
        
        WRITE(*,'(A,I0)') 'Test matrix dimension: n = ', n_analytical
        WRITE(*,'(A,I0)') 'Krylov size m = ', krylov_size
        WRITE(*,*) 'Expected eigenvalues (m largest):'
        DO k = 1, krylov_size
            WRITE(*,'(I3,2X,ES15.6,A,I0,A)') k, analytical_evals(k), '  (= ', n_analytical - k + 1, ')'
        END DO
        WRITE(*,*) ''
        
        ! Check 1: Are eigenvalues real? (Imaginary part should be ~0)
        max_imag = 0.0_rk
        DO k = 1, krylov_size
            max_imag = MAX(max_imag, ABS(AIMAG(eigenvalues(k))))
        END DO
        
        all_real = (max_imag < 1.0E-6_rk)
        WRITE(*,'(A,ES12.4)') 'Maximum imaginary part: ', max_imag
        IF (all_real) THEN
            WRITE(*,*) '✓ PASS: All eigenvalues are real (Im(λ) < 1E-6)'
        ELSE
            WRITE(*,*) '✗ FAIL: Some eigenvalues have significant imaginary parts'
        END IF
        WRITE(*,*) ''
        
        ! Check 2: Compare computed vs analytical eigenvalues
        ! For diagonal matrix, Arnoldi should capture the m largest exactly
        WRITE(*,*) 'Comparing computed vs analytical eigenvalues:'
        WRITE(*,*) '  #   Computed         Analytical       Error         % Error'
        WRITE(*,*) '---  --------------   --------------   ----------   -----------'
        
        max_error = 0.0_rk
        max_rel_error = 0.0_rk
        
        DO k = 1, krylov_size
            computed_val = REAL(eigenvalues(k))
            error_abs = ABS(computed_val - analytical_evals(k))
            rel_error = error_abs / analytical_evals(k) * 100.0_rk  ! Percentage
            
            max_error = MAX(max_error, error_abs)
            max_rel_error = MAX(max_rel_error, rel_error)
            
            WRITE(*,'(I3,2X,ES15.6,2X,ES15.6,2X,ES11.3,2X,F10.4,A)') k, &
                computed_val, analytical_evals(k), error_abs, rel_error, '%'
        END DO
        
        WRITE(*,*) '---  --------------   --------------   ----------   -----------'
        WRITE(*,'(A,ES12.4)') 'Maximum absolute error: ', max_error
        WRITE(*,'(A,F10.4,A)') 'Maximum relative error: ', max_rel_error, '%'
        WRITE(*,*) ''
        
        ! Check if errors are acceptable
        IF (max_error < 1.0E-6_rk) THEN
            WRITE(*,*) '✓ EXCELLENT: Eigenvalues match to machine precision!'
        ELSE IF (max_error < 1.0E-3_rk) THEN
            WRITE(*,*) '✓ GOOD: Eigenvalues match with < 0.1% error'
        ELSE IF (max_error < 0.1_rk) THEN
            WRITE(*,*) '✓ ACCEPTABLE: Eigenvalues match with < 10% error'
        ELSE
            WRITE(*,*) '✗ POOR: Eigenvalues have significant errors'
        END IF
        WRITE(*,*) ''
        
        ! Overall validation
        validation_passed = all_real .AND. (max_error < 1.0E-3_rk)
        
        IF (validation_passed) THEN
            WRITE(*,*) '==============================================='
            WRITE(*,*) '✓✓✓ VALIDATION PASSED ✓✓✓'
            WRITE(*,*) '==============================================='
            WRITE(*,*) 'Arnoldi computed eigenvalues are valid!'
            WRITE(*,*) '  - All eigenvalues are real'
            WRITE(*,*) '  - All eigenvalues within expected range'
        ELSE
            WRITE(*,*) '==============================================='
            WRITE(*,*) '✗✗✗ VALIDATION FAILED ✗✗✗'
            WRITE(*,*) '==============================================='
            IF (.NOT. all_real) THEN
                WRITE(*,*) '  - Eigenvalues have imaginary parts (should be real for test matrix)'
            END IF
            IF (max_error >= 0.1_rk) THEN
                WRITE(*,*) '  - Some eigenvalues outside expected range [-2, +2]'
            END IF
        END IF
        WRITE(*,*) ''
        
        DEALLOCATE(analytical_evals)
        
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
