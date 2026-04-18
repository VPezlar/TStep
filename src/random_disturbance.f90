! =============================================================================
! random_disturbance  --  generates a random, unit-norm perturbation vector
!                         scaled by a user-supplied magnitude.
!
! `initial_disturbance(N, mag, v, ierr)` returns an allocatable REAL vector v
! of length N with values drawn uniformly in [-1, 1), then L2-normalized and
! multiplied by `mag` (so ||v||_2 == mag). Used to seed the initial Krylov
! vector for Arnoldi and to produce small flowfield perturbations for the
! Frechet-derivative matvec.
!
! Relies on the intrinsic RANDOM_NUMBER generator (no explicit seeding yet
! so each run is reproducible with the compiler's default seed state).
! =============================================================================
MODULE random_disturbance
    USE accuracy
    USE variables
    USE setup, ONLY: INT_TO_STR
    USE error_handling

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: initial_disturbance


CONTAINS

    SUBROUTINE initial_disturbance(VECTOR_LENGTH, SCALING_CONSTANT, FINAL_VECTOR, ERROR_STATUS)
        ! Arguments
        INTEGER(ik), INTENT(IN)              :: VECTOR_LENGTH        ! The desired size of the output vector (N)
        REAL(rk), INTENT(IN)             :: SCALING_CONSTANT     ! The factor to multiply the normalized vector by
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: FINAL_VECTOR       ! The final real, normalized, and scaled vector
        INTEGER(ik), INTENT(OUT)             :: ERROR_STATUS         ! 0 for success, non-zero for failure

        ! Local Variables
        REAL(rk) :: NORM_VAL       ! The calculated L2-norm (magnitude) of the initial random vector
        INTEGER(ik) :: ALLOC_STAT     ! Status for allocation/deallocation

        ! --- Initialization ---
        ERROR_STATUS = 0

        ! Clean up existing array and validate input
        IF (ALLOCATED(FINAL_VECTOR)) THEN
            DEALLOCATE(FINAL_VECTOR, STAT=ALLOC_STAT)
        END IF

        IF (VECTOR_LENGTH <= 0) THEN
            ERROR_STATUS = ERR_DIST_INVALID_LENGTH
            CALL log_error(ERR_DIST_INVALID_LENGTH, 'Length specified: 0 or negative')
            RETURN
        END IF

        ! 1. Allocate the output vector
        ALLOCATE(FINAL_VECTOR(VECTOR_LENGTH), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_DIST_ALLOC
            CALL log_error(ERR_DIST_ALLOC, 'Allocation status: '//TRIM(ADJUSTL(INT_TO_STR(ALLOC_STAT))))
            RETURN
        END IF

        ! 2. Fill with uniform random REAL numbers in the range [0.0, 1.0)
        !    (These are the base values for normalization)
        CALL RANDOM_NUMBER(FINAL_VECTOR)

        ! 2b. Shift to [-1.0, 1.0) by scaling and translating
        FINAL_VECTOR = 2.0_rk * FINAL_VECTOR - 1.0_rk

        ! 3. Calculate the L2-Norm (Magnitude) of the current vector
        !    NORM2 is an intrinsic function in Fortran for the L2-norm (sqrt(sum(x_i^2)))
        NORM_VAL = NORM2(FINAL_VECTOR)

        ! Check if the norm is effectively zero (highly unlikely but good practice)
        IF (NORM_VAL < 1.0E-12_rk) THEN
            ERROR_STATUS = ERR_DIST_ZERO_NORM
            CALL log_error(ERR_DIST_ZERO_NORM, 'Computed norm is near zero')
            DEALLOCATE(FINAL_VECTOR)
            RETURN
        END IF

        ! 4. Normalize the vector (Divide by its norm)
        !    This makes the L2-norm of FINAL_VECTOR equal to 1.0
        FINAL_VECTOR = FINAL_VECTOR / NORM_VAL

        ! 5. Scale the normalized vector by the constant
        FINAL_VECTOR = FINAL_VECTOR * SCALING_CONSTANT

    END SUBROUTINE initial_disturbance

END MODULE random_disturbance
