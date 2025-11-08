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
        REAL(rk), DIMENSION(:), ALLOCATABLE :: computed_real_parts
        REAL(rk) :: max_error, avg_error, max_imag, error_val
        INTEGER(ik) :: k, n_check
        LOGICAL :: all_real, validation_passed
        
        WRITE(*,*) '==============================================='
        WRITE(*,*) 'EIGENVALUE VALIDATION'
        WRITE(*,*) '==============================================='
        WRITE(*,*) ''
        
        n_analytical = SIZE(v_normalized)
        pi = 4.0_rk * ATAN(1.0_rk)
        
        ! Compute analytical eigenvalues: λ_k = 2*cos(k*π/(n+1))
        ALLOCATE(analytical_evals(n_analytical))
        DO k = 1, n_analytical
            analytical_evals(k) = 2.0_rk * COS(REAL(k, rk) * pi / REAL(n_analytical + 1, rk))
        END DO
        
        WRITE(*,'(A,I0)') 'Test matrix dimension: n = ', n_analytical
        WRITE(*,*) 'Analytical eigenvalues (first 10):'
        DO k = 1, MIN(10, n_analytical)
            WRITE(*,'(I3,2X,ES15.6)') k, analytical_evals(k)
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
        
        ! Check 2: Compare real parts with analytical values
        ! Sort both arrays by magnitude for direct comparison
        n_check = MIN(krylov_size, 10)  ! Check first 10
        
        ALLOCATE(computed_real_parts(n_check))
        DO k = 1, n_check
            computed_real_parts(k) = REAL(eigenvalues(k))
        END DO
        
        ! Sort analytical eigenvalues by magnitude (descending)
        CALL sort_by_magnitude(analytical_evals, n_analytical)
        
        ! Sort computed eigenvalues by magnitude (descending)
        CALL sort_by_magnitude(computed_real_parts, n_check)
        
        WRITE(*,*) 'Comparing computed vs analytical eigenvalues:'
        WRITE(*,*) '(Both sorted by descending magnitude)'
        WRITE(*,*) '  #   Computed         Analytical       Error'
        WRITE(*,*) '---  --------------   --------------   ----------'
        
        max_error = 0.0_rk
        avg_error = 0.0_rk
        DO k = 1, n_check
            error_val = ABS(computed_real_parts(k) - analytical_evals(k))
            max_error = MAX(max_error, error_val)
            avg_error = avg_error + error_val
            
            WRITE(*,'(I3,2X,ES15.6,2X,ES15.6,2X,ES11.3)') k, &
                computed_real_parts(k), analytical_evals(k), error_val
        END DO
        avg_error = avg_error / REAL(n_check, rk)
        
        WRITE(*,*) '---  --------------   --------------   ----------'
        WRITE(*,'(A,ES12.4)') 'Maximum error: ', max_error
        WRITE(*,'(A,ES12.4)') 'Average error: ', avg_error
        WRITE(*,*) ''
        
        ! Overall validation
        validation_passed = all_real .AND. (max_error < 1.0E-3_rk)
        
        IF (validation_passed) THEN
            WRITE(*,*) '==============================================='
            WRITE(*,*) '✓✓✓ VALIDATION PASSED ✓✓✓'
            WRITE(*,*) '==============================================='
            WRITE(*,*) 'Computed eigenvalues match analytical values!'
        ELSE
            WRITE(*,*) '==============================================='
            WRITE(*,*) '✗✗✗ VALIDATION FAILED ✗✗✗'
            WRITE(*,*) '==============================================='
            WRITE(*,*) 'Computed eigenvalues differ from expected!'
            IF (.NOT. all_real) THEN
                WRITE(*,*) '  - Eigenvalues have imaginary parts (should be real)'
            END IF
            IF (max_error >= 1.0E-3_rk) THEN
                WRITE(*,*) '  - Errors exceed tolerance (max error > 1E-3)'
            END IF
        END IF
        WRITE(*,*) ''
        
        DEALLOCATE(analytical_evals, computed_real_parts)
        
    END SUBROUTINE validate_eigenvalues
    
    SUBROUTINE sort_by_magnitude(array, n)
        ! Simple bubble sort by descending magnitude
        REAL(rk), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER(ik), INTENT(IN) :: n
        INTEGER(ik) :: i, j
        REAL(rk) :: temp
        
        DO i = 1, n-1
            DO j = i+1, n
                IF (ABS(array(j)) > ABS(array(i))) THEN
                    temp = array(i)
                    array(i) = array(j)
                    array(j) = temp
                END IF
            END DO
        END DO
    END SUBROUTINE sort_by_magnitude
    
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
