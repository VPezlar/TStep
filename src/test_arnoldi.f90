PROGRAM test_arnoldi
    USE accuracy
    USE error_handling
    USE Arnoldi
    
    IMPLICIT NONE
    
    INTEGER(ik), PARAMETER :: n = 10   ! System dimension
    INTEGER(ik), PARAMETER :: m = 5    ! Krylov subspace size
    
    REAL(rk), DIMENSION(n) :: v_init
    COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: eigenvalues
    COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: eigenvectors
    INTEGER(ik) :: i
    INTEGER(ik) :: error_status
    CHARACTER(len=20) :: frechet_order
    REAL(rk) :: eps_0, TTime
    
    WRITE(*,*) '==============================================='
    WRITE(*,*) 'ARNOLDI MODULE - EIGENVALUE COMPUTATION TEST'
    WRITE(*,*) '==============================================='
    WRITE(*,*) ''
    
    ! Set parameters
    frechet_order = 'First'
    eps_0 = 1.0E-6_rk
    TTime = 1.0_rk
    
    ! Initialize random starting vector
    WRITE(*,'(A,I0)') 'System dimension n = ', n
    WRITE(*,'(A,I0)') 'Krylov size m = ', m
    WRITE(*,*) ''
    WRITE(*,*) 'Generating random initial disturbance vector...'
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(v_init)
    
    ! Normalize the vector (simulating what initial_disturbance module does)
    v_init = v_init / NORM2(v_init)
    WRITE(*,'(A,ES15.6)') 'Initial vector normalized: ||v_init|| = ', NORM2(v_init)
    
    ! Call Arnoldi eigenvalue routine
    WRITE(*,*) ''
    WRITE(*,*) 'Running Arnoldi iteration with eigenvalue computation...'
    WRITE(*,*) 'Using skip_normalization = .TRUE. (vector already normalized)'
    WRITE(*,*) ''
    CALL arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                            eigenvalues, eigenvectors, error_status, &
                            skip_normalization=.TRUE.)
    
    IF (error_status /= 0) THEN
        CALL log_error(error_status)
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
    DO i = 1, m
        WRITE(*,'(I3,2X,ES15.6,2X,ES15.6,2X,ES15.6)') i, &
            REAL(eigenvalues(i)), AIMAG(eigenvalues(i)), ABS(eigenvalues(i))
    END DO
    WRITE(*,*) '-----------------------------------------------'
    
    WRITE(*,*) ''
    WRITE(*,*) 'First Ritz vector (first 5 components):'
    DO i = 1, MIN(5, n)
        WRITE(*,'(I3,2X,ES15.6,2X,ES15.6)') i, &
            REAL(eigenvectors(i,1)), AIMAG(eigenvectors(i,1))
    END DO
    
    ! Cleanup
    IF (ALLOCATED(eigenvalues)) DEALLOCATE(eigenvalues)
    IF (ALLOCATED(eigenvectors)) DEALLOCATE(eigenvectors)
    
    WRITE(*,*) ''
    WRITE(*,*) '==============================================='
    WRITE(*,*) 'TEST COMPLETED SUCCESSFULLY'
    WRITE(*,*) '==============================================='
    
END PROGRAM test_arnoldi
