! =============================================================================
! frechet_stencil  --  computes stencil nodes and weights for the Frechet
!                      derivative finite-difference approximation.
!
! DF(q0)[v] ~= ( sum_i  w_i * F(q0 + alpha_i * eps_0 * v) ) / eps_0
!
! Nodes alpha_i and weights w_i solve the Vandermonde system:
!   V * w = e_1,   V(j,i) = alpha_i^(j-1),  e_1 has 1 in position 2.
!
! Supported orders:
!   order = 1     : two-point biased stencil  (N=2, nodes=[0,1])
!   order even    : symmetric stencil, N=order, nodes = {-m,...,-1,1,...,m}
!                   where m = order/2
!
! Public: Frechet_weights(order, nodes, weights, ierr)
! =============================================================================
MODULE frechet_stencil

    USE accuracy
    USE error_handling
    USE setup, ONLY: INT_TO_STR

    IMPLICIT NONE

    INTERFACE
        SUBROUTINE DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            IMPORT :: ik, rk
            INTEGER(ik), INTENT(IN)    :: N, NRHS, LDA, LDB
            INTEGER(ik), INTENT(OUT)   :: INFO
            INTEGER(ik), INTENT(OUT)   :: IPIV(*)
            REAL(rk),    INTENT(INOUT) :: A(LDA, *)
            REAL(rk),    INTENT(INOUT) :: B(LDB, *)
        END SUBROUTINE DGESV
    END INTERFACE

    PRIVATE
    PUBLIC :: Frechet_weights

CONTAINS

    SUBROUTINE Frechet_weights(order, nodes, weights, ierr)
        ! Arguments
        INTEGER(ik),                         INTENT(IN)  :: order
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nodes
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: weights
        INTEGER(ik),                         INTENT(OUT) :: ierr

        ! Local
        INTEGER(ik) :: N, icol, ideg, k, i
        INTEGER(ik) :: alloc_stat, info_lapack
        REAL(rk),    DIMENSION(:,:), ALLOCATABLE :: V
        INTEGER(ik), DIMENSION(:),   ALLOCATABLE :: ipiv

        ierr = 0

        ! --- Validate order and set stencil size N ---
        IF (order == 1) THEN
            N = 2
        ELSE IF (order >= 2 .AND. MOD(order, 2) == 0) THEN
            N = order
        ELSE
            ierr = ERR_FRECHET_INVALID_ORDER
            CALL log_error(ERR_FRECHET_INVALID_ORDER, &
                'Order must be 1 or a positive even integer, got: ' // &
                TRIM(ADJUSTL(INT_TO_STR(order))))
            RETURN
        END IF

        ! --- Allocate ---
        ALLOCATE(nodes(N), weights(N), V(N, N), ipiv(N), STAT=alloc_stat)
        IF (alloc_stat /= 0) THEN
            ierr = ERR_FRECHET_ALLOC
            CALL log_error(ERR_FRECHET_ALLOC, &
                'Allocation failed for stencil size N=' // TRIM(ADJUSTL(INT_TO_STR(N))))
            RETURN
        END IF

        ! --- Build nodes ---
        IF (order == 1) THEN
            nodes(1) = 0.0_rk
            nodes(2) = 1.0_rk
        ELSE
            k = 1
            DO i = -(order / 2), order / 2
                IF (i /= 0) THEN
                    nodes(k) = REAL(i, KIND=rk)
                    k = k + 1
                END IF
            END DO
        END IF

        ! --- Build Vandermonde matrix: V(ideg+1, icol) = nodes(icol)^ideg ---
        DO ideg = 0, N-1
            DO icol = 1, N
                V(ideg+1, icol) = nodes(icol)**ideg
            END DO
        END DO

        ! --- RHS: first-derivative condition (degree-1 slot = 1, rest = 0) ---
        weights    = 0.0_rk
        weights(2) = 1.0_rk

        ! --- Solve V * weights = RHS via LAPACK DGESV ---
        CALL DGESV(N, 1_ik, V, N, ipiv, weights, N, info_lapack)

        IF (info_lapack /= 0) THEN
            ierr = ERR_FRECHET_SOLVE
            CALL log_error(ERR_FRECHET_SOLVE, &
                'DGESV INFO=' // TRIM(ADJUSTL(INT_TO_STR(info_lapack))) // &
                ' order=' // TRIM(ADJUSTL(INT_TO_STR(order))))
            DEALLOCATE(nodes, weights, V, ipiv)
            RETURN
        END IF

        DEALLOCATE(V, ipiv)

    END SUBROUTINE Frechet_weights

END MODULE frechet_stencil
