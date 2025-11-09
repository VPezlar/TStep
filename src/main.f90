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
    
    ! Timing variables
    INTEGER :: clock_start, clock_end, clock_rate
    REAL(rk) :: elapsed_time

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
    ! Test matrix: n=1000, eigenvalues = [1, 2, 3, ..., 1000]
    ! Construction: A = eigvecs @ diag(eigvals) @ inv(eigvecs)
    ! Expected: Arnoldi recovers largest 100 eigenvalues: 1000, 999, ..., 901
    ! ===================================================================
    
    WRITE(*,*)
    WRITE(*,*) '======================================================='
    WRITE(*,*) 'ARNOLDI TEST - REFERENCE IMPLEMENTATION REPLICATION'
    WRITE(*,*) '======================================================='
    WRITE(*,*)
    WRITE(*,*) 'Test problem for scalability testing:'
    WRITE(*,*) '  - Dimension: n = 1000'
    WRITE(*,*) '  - Krylov size: m = 200 (large subspace for convergence)'
    WRITE(*,*) '  - Eigenvalues: 1, 2, 3, ..., 1000'
    WRITE(*,*) '  - Construction: A = Q @ diag(eigvals) @ Q^T (symmetric)'
    WRITE(*,*) '  - Q = orthogonal matrix from QR decomposition'
    WRITE(*,*)
    WRITE(*,*) 'Expected: Top 5 should be close to 1000, 999, 998, 997, 996'
    WRITE(*,*)
    
    ! Set number of threads for BLAS/LAPACK operations
    CALL set_blas_threads(num_threads)
    
    ! HARDCODED test size for scalability testing
    ALLOCATE(v_normalized(1000), STAT=error_status)
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
    
    ! Start timing
    CALL SYSTEM_CLOCK(clock_start, clock_rate)
    
    ! Call Arnoldi with n=1000, m=200 (larger subspace for better convergence)
    CALL arnoldi_eigenvalues(v_normalized, 200, frechet_order, eps_0, TTime, &
                            eigenvalues, eigenvectors, error_status, &
                            skip_normalization=.TRUE., sort_by='magnitude')
    
    ! End timing
    CALL SYSTEM_CLOCK(clock_end)
    elapsed_time = REAL(clock_end - clock_start, rk) / REAL(clock_rate, rk)
    
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
    WRITE(*,'(A,F12.3,A)') 'Arnoldi computation time: ', elapsed_time, ' seconds'
    WRITE(*,*)
    WRITE(*,*) 'Top 5 Ritz eigenvalues (sorted by magnitude):'
    WRITE(*,*) '-------------------------------------------------------'
    WRITE(*,*) '  #    Real Part      Imag Part      |λ|         Expected'
    WRITE(*,*) '-------------------------------------------------------'
    DO i = 1, MIN(5, SIZE(eigenvalues))
        WRITE(*,'(I3,2X,F12.6,2X,F12.6,2X,F12.6,2X,I5)') i, &
            REAL(eigenvalues(i)), AIMAG(eigenvalues(i)), ABS(eigenvalues(i)), 1001 - i
    END DO
    WRITE(*,*) '-------------------------------------------------------'
    WRITE(*,'(A,I0,A)') '(Computed ', SIZE(eigenvalues), ' eigenvalues total)'
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

    SUBROUTINE set_blas_threads(nthreads)
        ! Displays thread configuration for BLAS/LAPACK
        ! Note: Set environment variables BEFORE running the program:
        !   export OMP_NUM_THREADS=4
        !   export OPENBLAS_NUM_THREADS=4
        INTEGER(ik), INTENT(IN) :: nthreads
        
        IF (nthreads <= 0) THEN
            WRITE(*,*) 'BLAS threading: AUTO mode (using all available cores)'
            WRITE(*,*) 'To limit cores, set OMP_NUM_THREADS before running'
        ELSE
            WRITE(*,'(A,I0,A)') 'BLAS threading: Requesting ', nthreads, ' thread(s)'
            WRITE(*,*) 'Set these before running for OpenBLAS/MKL:'
            WRITE(*,'(A,I0)') '  export OMP_NUM_THREADS=', nthreads
            WRITE(*,'(A,I0)') '  export OPENBLAS_NUM_THREADS=', nthreads
        END IF
        WRITE(*,*)
        
    END SUBROUTINE set_blas_threads

    SUBROUTINE validate_eigenvalues()
        ! Validates test matrix eigenvalues: n=1000, m=200
        ! Matrix A has eigenvalues 1, 2, 3, ..., 1000
        ! With large Krylov subspace, top 5 should be very close to 1000, 999, 998, 997, 996
        REAL(rk) :: max_imag, max_eval, rel_error
        REAL(rk) :: computed_val, expected_val, max_rel_error
        INTEGER(ik) :: k, m_size
        LOGICAL :: all_real, top5_correct
        
        m_size = SIZE(eigenvalues)
        
        WRITE(*,*) '======================================================='
        WRITE(*,*) 'VALIDATION'
        WRITE(*,*) '======================================================='
        WRITE(*,*)
        
        ! Check 1: Are eigenvalues real?
        max_imag = 0.0_rk
        DO k = 1, m_size
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
        
        ! Check 2: Are top 5 eigenvalues correct?
        ! With m=200, top 5 should be very accurate (< 0.1% error)
        max_rel_error = 0.0_rk
        max_eval = 0.0_rk
        
        WRITE(*,*) 'Top 5 eigenvalues vs expected:'
        WRITE(*,*) '  #    Computed      Expected     Rel Error'
        WRITE(*,*) '-----------------------------------------------'
        DO k = 1, MIN(5, m_size)
            computed_val = ABS(eigenvalues(k))
            expected_val = REAL(1001 - k, rk)  ! 1000, 999, 998, 997, 996
            rel_error = ABS(computed_val - expected_val) / expected_val * 100.0_rk
            max_rel_error = MAX(max_rel_error, rel_error)
            max_eval = MAX(max_eval, computed_val)
            WRITE(*,'(I3,2X,F12.2,2X,F12.2,2X,F10.6,A)') k, computed_val, expected_val, rel_error, '%'
        END DO
        WRITE(*,*) '-----------------------------------------------'
        WRITE(*,*)
        
        ! Success if top 5 have < 1% error
        top5_correct = (max_rel_error < 1.0_rk)
        
        WRITE(*,'(A,F10.6,A)') 'Maximum relative error (top 5): ', max_rel_error, '%'
        IF (max_rel_error < 0.01_rk) THEN
            WRITE(*,*) 'EXCELLENT: Top 5 eigenvalues < 0.01% error'
        ELSE IF (max_rel_error < 0.1_rk) THEN
            WRITE(*,*) 'VERY GOOD: Top 5 eigenvalues < 0.1% error'
        ELSE IF (max_rel_error < 1.0_rk) THEN
            WRITE(*,*) 'GOOD: Top 5 eigenvalues < 1% error'
        ELSE
            WRITE(*,*) 'FAIL: Top 5 eigenvalues have > 1% error'
        END IF
        WRITE(*,*)
        
        ! Final verdict
        IF (all_real .AND. top5_correct) THEN
            WRITE(*,*) '======================================================='
            WRITE(*,*) 'SUCCESS: ARNOLDI WORKING CORRECTLY!'
            WRITE(*,*) '======================================================='
            WRITE(*,*)
            WRITE(*,*) 'Large Krylov subspace (m=200) successfully captures'
            WRITE(*,*) 'the dominant eigenvalues with high accuracy.'
            WRITE(*,*)
            WRITE(*,'(A,F6.3,A)') 'Baseline timing: ', elapsed_time, ' seconds'
            WRITE(*,*)
            WRITE(*,*) 'Ready for multi-threaded speedup testing!'
            WRITE(*,*) 'Run with different OMP_NUM_THREADS: 1, 2, 4, 8'
        ELSE
            WRITE(*,*) '======================================================='
            WRITE(*,*) 'FAILURE: IMPLEMENTATION HAS ISSUES'
            WRITE(*,*) '======================================================='
            IF (.NOT. all_real) THEN
                WRITE(*,*) '  Problem: Eigenvalues have imaginary components'
            END IF
            IF (.NOT. top5_correct) THEN
                WRITE(*,*) '  Problem: Top 5 eigenvalues not accurate enough'
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
