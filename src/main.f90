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
    WRITE(*,*) '  - Krylov size: m = 100'
    WRITE(*,*) '  - Eigenvalues: 1, 2, 3, ..., 1000'
    WRITE(*,*) '  - Construction: A = eigvecs @ diag(eigvals) @ inv(eigvecs)'
    WRITE(*,*) '  - Random eigenvector matrix (dense)'
    WRITE(*,*)
    WRITE(*,*) 'Expected: Largest 100 eigenvalues (1000, 999, ..., 901)'
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
    
    ! Call Arnoldi with n=1000, m=100
    CALL arnoldi_eigenvalues(v_normalized, 100, frechet_order, eps_0, TTime, &
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
    WRITE(*,*) 'Ritz eigenvalues (sorted by magnitude, showing first 10):'
    WRITE(*,*) '-------------------------------------------------------'
    WRITE(*,*) '  #    Real Part      Imag Part      |λ|'
    WRITE(*,*) '-------------------------------------------------------'
    DO i = 1, MIN(10, SIZE(eigenvalues))
        WRITE(*,'(I3,2X,F12.6,2X,F12.6,2X,F12.6)') i, &
            REAL(eigenvalues(i)), AIMAG(eigenvalues(i)), ABS(eigenvalues(i))
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
        ! Validates test matrix eigenvalues: n=1000, m=100
        ! Matrix A has eigenvalues 1, 2, 3, ..., 1000
        ! Arnoldi should find the 100 LARGEST (close to 1000)
        REAL(rk) :: max_imag, min_eval, max_eval, mean_eval
        REAL(rk) :: computed_val
        INTEGER(ik) :: k, m_size, num_in_range
        LOGICAL :: all_real, in_correct_range
        
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
        
        ! Check 2: Are eigenvalues in the correct range?
        ! Expected: largest 100 eigenvalues should be in range [850, 1000]
        min_eval = 1.0E10_rk
        max_eval = 0.0_rk
        mean_eval = 0.0_rk
        num_in_range = 0
        
        DO k = 1, m_size
            computed_val = ABS(eigenvalues(k))
            min_eval = MIN(min_eval, computed_val)
            max_eval = MAX(max_eval, computed_val)
            mean_eval = mean_eval + computed_val
            IF (computed_val >= 850.0_rk .AND. computed_val <= 1000.0_rk) THEN
                num_in_range = num_in_range + 1
            END IF
        END DO
        mean_eval = mean_eval / REAL(m_size, rk)
        
        WRITE(*,'(A,F10.2)') 'Largest eigenvalue:  ', max_eval
        WRITE(*,'(A,F10.2)') 'Smallest eigenvalue: ', min_eval
        WRITE(*,'(A,F10.2)') 'Mean eigenvalue:     ', mean_eval
        WRITE(*,'(A,I0,A,I0)') 'Eigenvalues in range [850,1000]: ', num_in_range, '/', m_size
        WRITE(*,*)
        
        ! Success if: max > 990 and mean > 900 and most are in range
        in_correct_range = (max_eval > 990.0_rk) .AND. &
                          (mean_eval > 900.0_rk) .AND. &
                          (num_in_range >= INT(0.8_rk * m_size))
        
        IF (in_correct_range) THEN
            WRITE(*,*) 'PASS: Eigenvalues are in the expected range (largest 100)'
        ELSE
            WRITE(*,*) 'FAIL: Eigenvalues not in expected range'
            WRITE(*,*) '      Expected: max>990, mean>900, 80%+ in [850,1000]'
        END IF
        WRITE(*,*)
        
        ! Final verdict
        IF (all_real .AND. in_correct_range) THEN
            WRITE(*,*) '======================================================='
            WRITE(*,*) 'SUCCESS: ARNOLDI IMPLEMENTATION IS WORKING!'
            WRITE(*,*) '======================================================='
            WRITE(*,*)
            WRITE(*,*) 'Arnoldi successfully finds the largest eigenvalues.'
            WRITE(*,'(A,F6.3,A)') 'Ready for scalability testing (time: ', elapsed_time, 's)'
        ELSE
            WRITE(*,*) '======================================================='
            WRITE(*,*) 'FAILURE: IMPLEMENTATION HAS ISSUES'
            WRITE(*,*) '======================================================='
            IF (.NOT. all_real) THEN
                WRITE(*,*) '  Problem: Eigenvalues have imaginary components'
            END IF
            IF (.NOT. in_correct_range) THEN
                WRITE(*,*) '  Problem: Eigenvalues not in expected range'
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
