! =============================================================================
! test_state_vector  --  standalone unit test for pack_state / unpack_state.
!
! Verifies:
!   1. state_length(Nmesh) returns NVARS*Nmesh.
!   2. Pack followed by unpack recovers the original fields bit-for-bit.
!   3. Pack of known-pattern fields places each field in the expected
!      block of the state vector (layout contract).
!   4. Shape-mismatch errors are detected.
!
! Build (adjust paths to your tree):
!   gfortran -c ../src/accuracy.f90 -o ../obj/accuracy.o -J ../mod/
!   gfortran -c ../src/error_handling.f90 -o ../obj/error_handling.o -J ../mod/
!   gfortran -c ../src/variables.f90 -o ../obj/variables.o -J ../mod/
!   gfortran -c ../src/setup.f90 -o ../obj/setup.o -J ../mod/
!   gfortran -c ../src/state_vector.f90 -o ../obj/state_vector.o -J ../mod/
!   gfortran -c ../src/test_state_vector.f90 -o ../obj/test_state_vector.o -J ../mod/
!   gfortran ../obj/*.o -o test_state_vector
!
! Run:   ./test_state_vector
! Expect: ALL TESTS PASSED
! =============================================================================
PROGRAM test_state_vector
    USE accuracy
    USE state_vector
    IMPLICIT NONE

    INTEGER(ik), PARAMETER :: Nmesh = 5_ik
    INTEGER(ik)            :: N

    REAL(rk), DIMENSION(Nmesh) :: rho_in,  p_in,  T_in,  U_in,  V_in,  W_in
    REAL(rk), DIMENSION(Nmesh) :: rho_out, p_out, T_out, U_out, V_out, W_out
    REAL(rk), DIMENSION(:), ALLOCATABLE :: state_vec

    INTEGER(ik) :: ierr, i
    LOGICAL     :: all_passed
    REAL(rk), PARAMETER :: TOL = 0.0_rk   ! bit-for-bit identity expected

    all_passed = .TRUE.

    WRITE(*,'(A)') '==============================================='
    WRITE(*,'(A)') 'state_vector unit test'
    WRITE(*,'(A)') '==============================================='

    ! --- Test 1: state_length -------------------------------------------------
    N = state_length(Nmesh)
    IF (N == NVARS * Nmesh) THEN
        WRITE(*,'(A,I0)') '  [PASS] state_length(Nmesh) = NVARS*Nmesh = ', N
    ELSE
        WRITE(*,'(A,I0,A,I0)') '  [FAIL] state_length(Nmesh)=', N, &
                               ' expected ', NVARS*Nmesh
        all_passed = .FALSE.
    END IF

    ! --- Test 2: round-trip on known numeric fields ---------------------------
    !   Use a distinctive pattern per variable so mis-ordering shows up loud.
    DO i = 1, Nmesh
        rho_in(i) = 100.0_rk + REAL(i, rk)      ! 101, 102, 103, 104, 105
        p_in  (i) = 200.0_rk + REAL(i, rk)      ! 201, 202, 203, 204, 205
        T_in  (i) = 300.0_rk + REAL(i, rk)
        U_in  (i) = 400.0_rk + REAL(i, rk)
        V_in  (i) = 500.0_rk + REAL(i, rk)
        W_in  (i) = 600.0_rk + REAL(i, rk)
    END DO

    ALLOCATE(state_vec(N))

    CALL pack_state(rho_in, p_in, T_in, U_in, V_in, W_in, state_vec, ierr)
    IF (ierr /= 0) THEN
        WRITE(*,'(A,I0)') '  [FAIL] pack_state returned ierr=', ierr
        all_passed = .FALSE.
    END IF

    CALL unpack_state(state_vec, rho_out, p_out, T_out, U_out, V_out, W_out, ierr)
    IF (ierr /= 0) THEN
        WRITE(*,'(A,I0)') '  [FAIL] unpack_state returned ierr=', ierr
        all_passed = .FALSE.
    END IF

    IF ( ALL(rho_out == rho_in) .AND. ALL(p_out == p_in) .AND. &
         ALL(T_out   == T_in)   .AND. ALL(U_out == U_in) .AND. &
         ALL(V_out   == V_in)   .AND. ALL(W_out == W_in) ) THEN
        WRITE(*,'(A)') '  [PASS] pack then unpack reproduces all fields bit-for-bit'
    ELSE
        WRITE(*,'(A)') '  [FAIL] round-trip mismatch'
        WRITE(*,'(A,6(1X,F8.2))') '     rho_in :', rho_in
        WRITE(*,'(A,6(1X,F8.2))') '     rho_out:', rho_out
        all_passed = .FALSE.
    END IF

    ! --- Test 3: layout contract ----------------------------------------------
    !   rho should be in slots 1..Nmesh, p in Nmesh+1..2*Nmesh, etc.
    IF ( state_vec(1)             == 101.0_rk .AND. &
         state_vec(Nmesh + 1)     == 201.0_rk .AND. &
         state_vec(2*Nmesh + 1)   == 301.0_rk .AND. &
         state_vec(3*Nmesh + 1)   == 401.0_rk .AND. &
         state_vec(4*Nmesh + 1)   == 501.0_rk .AND. &
         state_vec(5*Nmesh + 1)   == 601.0_rk ) THEN
        WRITE(*,'(A)') '  [PASS] block layout rho,p,T,U,V,W in that order'
    ELSE
        WRITE(*,'(A)') '  [FAIL] state vector layout mismatch'
        WRITE(*,'(A)') '     first entry of each block:'
        WRITE(*,'(A,F8.2,A)') '       slot 1            =', state_vec(1),           '  (expected 101.00 for rho)'
        WRITE(*,'(A,F8.2,A)') '       slot Nmesh+1      =', state_vec(Nmesh+1),     '  (expected 201.00 for p)'
        WRITE(*,'(A,F8.2,A)') '       slot 2*Nmesh+1    =', state_vec(2*Nmesh+1),   '  (expected 301.00 for T)'
        WRITE(*,'(A,F8.2,A)') '       slot 3*Nmesh+1    =', state_vec(3*Nmesh+1),   '  (expected 401.00 for U)'
        WRITE(*,'(A,F8.2,A)') '       slot 4*Nmesh+1    =', state_vec(4*Nmesh+1),   '  (expected 501.00 for V)'
        WRITE(*,'(A,F8.2,A)') '       slot 5*Nmesh+1    =', state_vec(5*Nmesh+1),   '  (expected 601.00 for W)'
        all_passed = .FALSE.
    END IF

    ! --- Test 4: shape-mismatch detection -------------------------------------
    !   Call pack with a state_vec that is one element too short; expect ierr/=0.
    DEALLOCATE(state_vec)
    ALLOCATE(state_vec(N - 1))
    CALL pack_state(rho_in, p_in, T_in, U_in, V_in, W_in, state_vec, ierr)
    IF (ierr /= 0) THEN
        WRITE(*,'(A,I0,A)') '  [PASS] wrong-size state_vec correctly rejected (ierr=', ierr, ')'
    ELSE
        WRITE(*,'(A)') '  [FAIL] wrong-size state_vec NOT detected'
        all_passed = .FALSE.
    END IF
    DEALLOCATE(state_vec)

    ! --- Summary --------------------------------------------------------------
    WRITE(*,'(A)') '==============================================='
    IF (all_passed) THEN
        WRITE(*,'(A)') 'ALL TESTS PASSED'
        STOP 0
    ELSE
        WRITE(*,'(A)') 'SOME TESTS FAILED'
        STOP 1
    END IF

END PROGRAM test_state_vector