! =============================================================================
! state_vector  --  pack and unpack between physical-field representation
!                   (six per-cell arrays) and the flat state vector used by
!                   Arnoldi.
!
! The rest of TStep speaks in six 1-D arrays of length Nmesh:
!   rho, p, T, U, V, W
! Arnoldi needs ONE contiguous vector of length NVARS * Nmesh = 6 * Nmesh.
!
! Layout (block, not interleaved):
!   state_vec(          1 :   Nmesh) = rho(1:Nmesh)
!   state_vec(  Nmesh + 1 : 2*Nmesh) = p  (1:Nmesh)
!   state_vec(2*Nmesh + 1 : 3*Nmesh) = T  (1:Nmesh)
!   state_vec(3*Nmesh + 1 : 4*Nmesh) = U  (1:Nmesh)
!   state_vec(4*Nmesh + 1 : 5*Nmesh) = V  (1:Nmesh)
!   state_vec(5*Nmesh + 1 : 6*Nmesh) = W  (1:Nmesh)
!
! Keep this layout consistent EVERYWHERE a state vector is packed or
! unpacked. Mismatches silently corrupt the spectrum.
!
! Public:
!   NVARS                    -- parameter, number of per-cell variables (6)
!   state_length(Nmesh)      -- returns NVARS * Nmesh
!   pack_state  (rho, p, T, U, V, W, state_vec, ierr)
!   unpack_state(state_vec, rho, p, T, U, V, W, ierr)
! =============================================================================
MODULE state_vector

    USE accuracy
    USE error_handling
    USE setup, ONLY: INT_TO_STR

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: NVARS, state_length, pack_state, unpack_state

    ! Number of per-cell variables packed into the state vector.
    ! Order (canonical):  1=rho, 2=p, 3=T, 4=U, 5=V, 6=W
    INTEGER(ik), PARAMETER :: NVARS = 6_ik

CONTAINS

    ! -------------------------------------------------------------------------
    ! state_length  --  convenience helper for total state vector length.
    ! -------------------------------------------------------------------------
    PURE FUNCTION state_length(Nmesh) RESULT(N)
        INTEGER(ik), INTENT(IN) :: Nmesh
        INTEGER(ik)             :: N
        N = NVARS * Nmesh
    END FUNCTION state_length


    ! -------------------------------------------------------------------------
    ! pack_state  --  copy six per-cell arrays into one flat state vector.
    !
    ! All six input arrays must have the same length (Nmesh). The output
    ! vector must be already allocated to length NVARS * Nmesh. The routine
    ! does NOT allocate; the caller owns the memory. This lets pack_state
    ! be called in a tight loop (e.g. inside the Fréchet stencil) without
    ! repeated alloc/dealloc traffic.
    ! -------------------------------------------------------------------------
    SUBROUTINE pack_state(rho, p, T, U, V, W, state_vec, ierr)
        REAL(rk), DIMENSION(:), INTENT(IN)  :: rho, p, T, U, V, W
        REAL(rk), DIMENSION(:), INTENT(OUT) :: state_vec
        INTEGER(ik),            INTENT(OUT) :: ierr

        INTEGER(ik) :: Nmesh, expected_len

        ierr = 0

        Nmesh = SIZE(rho)

        ! --- Shape checks (cheap but catch layout bugs immediately) ---
        IF (SIZE(p) /= Nmesh .OR. SIZE(T) /= Nmesh .OR. &
            SIZE(U) /= Nmesh .OR. SIZE(V) /= Nmesh .OR. SIZE(W) /= Nmesh) THEN
            ierr = ERR_STATE_SHAPE_MISMATCH
            CALL log_error(ERR_STATE_SHAPE_MISMATCH, &
                'pack: field arrays differ in length. ' // &
                'rho=' // TRIM(ADJUSTL(INT_TO_STR(Nmesh)))   // &
                ' p='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(p)))) // &
                ' T='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(T)))) // &
                ' U='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(U)))) // &
                ' V='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(V)))) // &
                ' W='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(W)))))
            RETURN
        END IF

        expected_len = NVARS * Nmesh
        IF (SIZE(state_vec) /= expected_len) THEN
            ierr = ERR_STATE_VECTOR_SIZE
            CALL log_error(ERR_STATE_VECTOR_SIZE, &
                'pack: state_vec length ' // &
                TRIM(ADJUSTL(INT_TO_STR(SIZE(state_vec)))) // &
                ' /= NVARS*Nmesh = ' // &
                TRIM(ADJUSTL(INT_TO_STR(expected_len))))
            RETURN
        END IF

        ! --- Block layout: six contiguous chunks of length Nmesh ---
        state_vec(           1 :   Nmesh) = rho
        state_vec(  Nmesh  + 1 : 2*Nmesh) = p
        state_vec(2*Nmesh  + 1 : 3*Nmesh) = T
        state_vec(3*Nmesh  + 1 : 4*Nmesh) = U
        state_vec(4*Nmesh  + 1 : 5*Nmesh) = V
        state_vec(5*Nmesh  + 1 : 6*Nmesh) = W

    END SUBROUTINE pack_state


    ! -------------------------------------------------------------------------
    ! unpack_state  --  copy one flat state vector back into six per-cell
    !                   arrays. Inverse of pack_state.
    !
    ! All six output arrays must be already allocated to the same length
    ! (Nmesh). The routine does NOT allocate; the caller owns the memory.
    ! -------------------------------------------------------------------------
    SUBROUTINE unpack_state(state_vec, rho, p, T, U, V, W, ierr)
        REAL(rk), DIMENSION(:), INTENT(IN)  :: state_vec
        REAL(rk), DIMENSION(:), INTENT(OUT) :: rho, p, T, U, V, W
        INTEGER(ik),            INTENT(OUT) :: ierr

        INTEGER(ik) :: Nmesh, expected_len

        ierr = 0

        Nmesh = SIZE(rho)

        ! --- Shape checks ---
        IF (SIZE(p) /= Nmesh .OR. SIZE(T) /= Nmesh .OR. &
            SIZE(U) /= Nmesh .OR. SIZE(V) /= Nmesh .OR. SIZE(W) /= Nmesh) THEN
            ierr = ERR_STATE_SHAPE_MISMATCH
            CALL log_error(ERR_STATE_SHAPE_MISMATCH, &
                'unpack: field arrays differ in length. ' // &
                'rho=' // TRIM(ADJUSTL(INT_TO_STR(Nmesh)))   // &
                ' p='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(p)))) // &
                ' T='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(T)))) // &
                ' U='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(U)))) // &
                ' V='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(V)))) // &
                ' W='   // TRIM(ADJUSTL(INT_TO_STR(SIZE(W)))))
            RETURN
        END IF

        expected_len = NVARS * Nmesh
        IF (SIZE(state_vec) /= expected_len) THEN
            ierr = ERR_STATE_VECTOR_SIZE
            CALL log_error(ERR_STATE_VECTOR_SIZE, &
                'unpack: state_vec length ' // &
                TRIM(ADJUSTL(INT_TO_STR(SIZE(state_vec)))) // &
                ' /= NVARS*Nmesh = ' // &
                TRIM(ADJUSTL(INT_TO_STR(expected_len))))
            RETURN
        END IF

        ! --- Reverse of pack_state (same layout) ---
        rho = state_vec(           1 :   Nmesh)
        p   = state_vec(  Nmesh  + 1 : 2*Nmesh)
        T   = state_vec(2*Nmesh  + 1 : 3*Nmesh)
        U   = state_vec(3*Nmesh  + 1 : 4*Nmesh)
        V   = state_vec(4*Nmesh  + 1 : 5*Nmesh)
        W   = state_vec(5*Nmesh  + 1 : 6*Nmesh)

    END SUBROUTINE unpack_state

END MODULE state_vector