PROGRAM test_arnoldi
    USE accuracy
    USE setup
    USE Arnoldi
    
    IMPLICIT NONE
    
    INTEGER(ik), PARAMETER :: n = 5    ! Matrix dimension
    INTEGER(ik), PARAMETER :: m = 3    ! Krylov subspace size
    
    REAL(rk), DIMENSION(n, n)   :: A
    REAL(rk), DIMENSION(n, m+1) :: v
    REAL(rk), DIMENSION(m+1, m) :: H
    INTEGER(ik) :: i, j
    
    WRITE(*,*) "Testing Arnoldi module integration..."
    
    ! Create a simple test matrix
    A = 0.0_rk
    DO i = 1, n
        A(i, i) = REAL(i, rk)
        IF (i < n) A(i, i+1) = 1.0_rk
    END DO
    
    ! Initialize first vector with random values
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(v(:,1))
    
    ! Call Arnoldi iteration
    WRITE(*,*) "Calling arnoldi_iter..."
    CALL arnoldi_iter(A, n, m, v, H)
    
    WRITE(*,*) "Arnoldi iteration completed successfully!"
    WRITE(*,*) "Hessenberg matrix H(1:3,1:3):"
    DO i = 1, 3
        WRITE(*,'(3ES15.6)') (H(i,j), j=1,3)
    END DO
    
END PROGRAM test_arnoldi
