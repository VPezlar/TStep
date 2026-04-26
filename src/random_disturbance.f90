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

    SUBROUTINE initial_disturbance(N_mesh, SCALING_CONSTANT, &
                                   rho0, p0, T0, &
                                   FINAL_VECTOR, ERROR_STATUS)
        ! Generates a thermodynamically-consistent random Krylov vector for the
        ! initial Arnoldi disturbance. The pressure block is constrained to satisfy
        ! the linearized perfect-gas EOS cell-by-cell:
        !
        !     p'(i) = (p0(i)/rho0(i)) * rho'(i)  +  (p0(i)/T0(i)) * T'(i)
        !
        ! This avoids the silent thermodynamic-projection that hePsiThermo applies
        ! to inconsistent input states (it would otherwise overwrite our chosen
        ! rho perturbation, wasting Krylov rank).
        !
        ! The gas constant R is implicitly captured per-cell via R(i) = p0/(rho0*T0)
        ! so this works regardless of dimensional or non-dimensional convention.

        INTEGER(ik),                    INTENT(IN)  :: N_mesh
        REAL(rk),                       INTENT(IN)  :: SCALING_CONSTANT
        REAL(rk),    DIMENSION(:),      INTENT(IN)  :: rho0, p0, T0
        REAL(rk),    DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: FINAL_VECTOR
        INTEGER(ik),                    INTENT(OUT) :: ERROR_STATUS

        INTEGER(ik) :: vlen, i, alloc_stat
        REAL(rk)    :: norm_val
        REAL(rk)    :: rho_prime, T_prime, p_prime
        INTEGER(ik) :: i_rho, i_p, i_t

        ERROR_STATUS = 0

        ! Validate
        IF (N_mesh <= 0) THEN
            ERROR_STATUS = ERR_DIST_INVALID_LENGTH
            CALL log_error(ERR_DIST_INVALID_LENGTH, 'N_mesh <= 0')
            RETURN
        END IF

        IF (SIZE(rho0) /= N_mesh .OR. SIZE(p0) /= N_mesh .OR. &
            SIZE(T0)   /= N_mesh) THEN
            ERROR_STATUS = ERR_DIST_INVALID_LENGTH
            CALL log_error(ERR_DIST_INVALID_LENGTH, &
                'Base-flow array size mismatch')
            RETURN
        END IF

        vlen = state_length(N_mesh)

        IF (ALLOCATED(FINAL_VECTOR)) DEALLOCATE(FINAL_VECTOR)
        ALLOCATE(FINAL_VECTOR(vlen), STAT=alloc_stat)
        IF (alloc_stat /= 0) THEN
            ERROR_STATUS = ERR_DIST_ALLOC
            CALL log_error(ERR_DIST_ALLOC, 'allocation failed')
            RETURN
        END IF

        ! Optional: explicit deterministic seed for cross-compiler reproducibility
        BLOCK
            INTEGER :: seed_size
            INTEGER, ALLOCATABLE :: seed(:)
            CALL RANDOM_SEED(size=seed_size)
            ALLOCATE(seed(seed_size))
            seed = 42
            CALL RANDOM_SEED(put=seed)
            DEALLOCATE(seed)
        END BLOCK

        ! 1. Generate uniform [-1, 1] noise for all 6 blocks
        CALL RANDOM_NUMBER(FINAL_VECTOR)
        FINAL_VECTOR = 2.0_rk * FINAL_VECTOR - 1.0_rk

        ! 2. Overwrite the p block with the thermodynamically consistent value
        DO i = 1, N_mesh
            i_rho = i
            i_p   = N_mesh + i
            i_t   = 2*N_mesh + i

            ! Sanity: avoid division-by-zero
            IF (rho0(i) < 1.0E-30_rk .OR. T0(i) < 1.0E-30_rk) THEN
                ERROR_STATUS = ERR_DIST_INVALID_LENGTH
                CALL log_error(ERR_DIST_INVALID_LENGTH, &
                    'rho0 or T0 vanishes -- bad base flow')
                DEALLOCATE(FINAL_VECTOR)
                RETURN
            END IF

            rho_prime = FINAL_VECTOR(i_rho)
            T_prime   = FINAL_VECTOR(i_t)

            ! Linearized EOS: p' = (p0/rho0) rho' + (p0/T0) T'
            p_prime = (p0(i) / rho0(i)) * rho_prime + &
                      (p0(i) / T0(i))   * T_prime

            FINAL_VECTOR(i_p) = p_prime
        END DO

        ! 3. Normalize to unit norm, then scale by mag
        norm_val = NORM2(FINAL_VECTOR)
        IF (norm_val < 1.0E-12_rk) THEN
            ERROR_STATUS = ERR_DIST_ZERO_NORM
            CALL log_error(ERR_DIST_ZERO_NORM, 'norm too small')
            DEALLOCATE(FINAL_VECTOR)
            RETURN
        END IF
        FINAL_VECTOR = FINAL_VECTOR / norm_val
        FINAL_VECTOR = FINAL_VECTOR * SCALING_CONSTANT

    END SUBROUTINE initial_disturbance
END MODULE random_disturbance