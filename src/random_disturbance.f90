MODULE random_disturbance
    USE accuracy
    USE variables
    USE setup

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: initial_disturbance


CONTAINS

    SUBROUTINE initial_disturbance(VECTOR_LENGTH, SCALING_CONSTANT, FINAL_VECTOR, ERROR_STATUS)
        ! Arguments
        INTEGER, INTENT(IN)              :: VECTOR_LENGTH        ! The desired size of the output vector (N)
        REAL(rk), INTENT(IN)             :: SCALING_CONSTANT     ! The factor to multiply the normalized vector by
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: FINAL_VECTOR       ! The final real, normalized, and scaled vector
        INTEGER, INTENT(OUT)             :: ERROR_STATUS         ! 0 for success, non-zero for failure

        ! Local Variables
        REAL(rk) :: NORM_VAL       ! The calculated L2-norm (magnitude) of the initial random vector
        INTEGER :: ALLOC_STAT     ! Status for allocation/deallocation

        ! --- Initialization ---
        ERROR_STATUS = 0

        ! Clean up existing array and validate input
        IF (ALLOCATED(FINAL_VECTOR)) THEN
            DEALLOCATE(FINAL_VECTOR, STAT=ALLOC_STAT)
        END IF

        IF (VECTOR_LENGTH <= 0) THEN
            WRITE(*,*) 'FATAL: Vector length must be positive.'
            ERROR_STATUS = 1
            RETURN
        END IF

        ! 1. Allocate the output vector
        ALLOCATE(FINAL_VECTOR(VECTOR_LENGTH), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            WRITE(*,*) 'FATAL: Failed to allocate memory for the final vector.'
            ERROR_STATUS = 2
            RETURN
        END IF

        ! 2. Fill with uniform random REAL numbers in the range [0.0, 1.0)
        !    (These are the base values for normalization)
        CALL RANDOM_NUMBER(FINAL_VECTOR)

        ! 3. Calculate the L2-Norm (Magnitude) of the current vector
        !    NORM2 is an intrinsic function in Fortran for the L2-norm (sqrt(sum(x_i^2)))
        NORM_VAL = NORM2(FINAL_VECTOR)

        ! Check if the norm is effectively zero (highly unlikely but good practice)
        IF (NORM_VAL < 1.0E-12_rk) THEN
            WRITE(*,*) 'WARNING: Norm is zero or near-zero. Cannot normalize.'
            DEALLOCATE(FINAL_VECTOR)
            ERROR_STATUS = 3
            RETURN
        END IF

        ! 4. Normalize the vector (Divide by its norm)
        !    This makes the L2-norm of FINAL_VECTOR equal to 1.0
        FINAL_VECTOR = FINAL_VECTOR / NORM_VAL

        ! 5. Scale the normalized vector by the constant
        FINAL_VECTOR = FINAL_VECTOR * SCALING_CONSTANT

    END SUBROUTINE initial_disturbance

END MODULE random_disturbance