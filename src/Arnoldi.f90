MODULE Arnoldi
    USE accuracy
    USE error_handling
    USE variables

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: arnoldi_eigenvalues

CONTAINS

    SUBROUTINE arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                                   eigenvalues, eigenvectors, ERROR_STATUS, &
                                   skip_normalization, sort_by)
        ! Computes Ritz eigenvalues and eigenvectors of the linearized flow
        ! operator using Arnoldi iteration.
        !
        ! The matrix-vector product w = A * v is supplied by the internal
        ! routine apply_linearized_operator, which is intended to approximate
        ! A via a Frechet-derivative finite difference using external CFD
        ! solver calls.
        !
        ! TODO: apply_linearized_operator is currently a STUB that returns
        !       ERR_ARNOLDI_NOT_IMPLEMENTED. arnoldi_eigenvalues will therefore
        !       fail on the first iteration until the real matvec is wired up.
        !
        ! If skip_normalization=.TRUE., the caller MUST ensure ||v_init|| = 1.

        ! Arguments
        REAL(rk), DIMENSION(:),              INTENT(IN)  :: v_init
        INTEGER(ik),                         INTENT(IN)  :: m
        CHARACTER(len=*),                    INTENT(IN)  :: frechet_order
        REAL(rk),                            INTENT(IN)  :: eps_0
        REAL(rk),                            INTENT(IN)  :: TTime
        COMPLEX(rk), DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: eigenvalues
        COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: eigenvectors
        INTEGER(ik),                         INTENT(OUT) :: ERROR_STATUS
        LOGICAL,          OPTIONAL,          INTENT(IN)  :: skip_normalization
        CHARACTER(len=*), OPTIONAL,          INTENT(IN)  :: sort_by

        ! Local Variables
        INTEGER(ik) :: n
        REAL(rk), DIMENSION(:,:), ALLOCATABLE :: V, H, H_m
        REAL(rk), DIMENSION(:),   ALLOCATABLE :: w, w_temp
        REAL(rk) :: norm_w, h_correction
        INTEGER(ik) :: i, j, k
        INTEGER(ik) :: ALLOC_STAT
        LOGICAL     :: do_normalize

        ! LAPACK and sorting workspaces
        COMPLEX(rk), DIMENSION(:),   ALLOCATABLE :: eval_work
        COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: evec_work
        REAL(rk),    DIMENSION(:),   ALLOCATABLE :: RWORK
        INTEGER(ik), DIMENSION(:),   ALLOCATABLE :: sort_idx
        REAL(rk),    DIMENSION(:),   ALLOCATABLE :: imag_parts
        COMPLEX(rk), DIMENSION(:),   ALLOCATABLE :: temp_evec

        ! --- Initialization ---
        ERROR_STATUS = 0
        n = SIZE(v_init)

        IF (n <= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_INVALID_DIM
            CALL log_error(ERR_ARNOLDI_INVALID_DIM, 'Initial vector has zero size')
            RETURN
        END IF

        IF (m <= 0 .OR. m > n) THEN
            ERROR_STATUS = ERR_ARNOLDI_INVALID_KRYLOV
            CALL log_error(ERR_ARNOLDI_INVALID_KRYLOV, 'Invalid Krylov size m')
            RETURN
        END IF

        ! --- Allocate arrays ---
        ALLOCATE(V(n, m+1), H(m+1, m), H_m(m, m), w(n), w_temp(n), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate Arnoldi arrays')
            RETURN
        END IF

        ALLOCATE(eval_work(m), evec_work(n, m), RWORK(2*m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate workspace arrays')
            DEALLOCATE(V, H, H_m, w, w_temp)
            RETURN
        END IF

        ! --- Step 1: Initialize first Krylov vector from v_init ---
        do_normalize = .TRUE.
        IF (PRESENT(skip_normalization)) THEN
            IF (skip_normalization) do_normalize = .FALSE.
        END IF

        V(:, 1) = v_init
        IF (do_normalize) THEN
            norm_w = NORM2(V(:, 1))
            IF (norm_w < 1.0E-14_rk) THEN
                ERROR_STATUS = ERR_ARNOLDI_ZERO_V1
                CALL log_error(ERR_ARNOLDI_ZERO_V1, 'Initial vector has zero norm')
                DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK)
                RETURN
            END IF
            V(:, 1) = V(:, 1) / norm_w
        END IF

        H = 0.0_rk

        ! --- Step 2: Arnoldi iteration ---
        DO j = 1, m
            ! w = A * v_j via the linearized-operator matvec (CFD-based)
            CALL apply_linearized_operator(V(:, j), w, frechet_order, eps_0, TTime, ERROR_STATUS)
            IF (ERROR_STATUS /= 0) THEN
                ! Matvec failed (or is still the stub) — abort cleanly.
                DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK)
                RETURN
            END IF

            ! Classical Gram-Schmidt (first pass)
            DO i = 1, j
                H(i, j) = DOT_PRODUCT(V(:, i), w)
                w = w - H(i, j) * V(:, i)
            END DO

            ! Reorthogonalization pass (critical for numerical stability with large m)
            DO i = 1, j
                h_correction = DOT_PRODUCT(V(:, i), w)
                H(i, j) = H(i, j) + h_correction
                w = w - h_correction * V(:, i)
            END DO

            norm_w = NORM2(w)
            H(j+1, j) = norm_w

            IF (norm_w > 1.0E-12_rk) THEN
                V(:, j+1) = w / norm_w
            ELSE
                WRITE(*,'(A,I0,A)') 'WARNING: Arnoldi breakdown at iteration ', j, &
                                    '. Krylov subspace exhausted.'
                V(:, j+1:m+1) = 0.0_rk
                H(j+2:m+1, :) = 0.0_rk
                EXIT
            END IF
        END DO

        ! --- Step 3: Extract m×m upper block of Hessenberg matrix ---
        H_m = H(1:m, 1:m)

        ! --- Step 4: Eigenvalues and eigenvectors of H_m (LAPACK DGEEV) ---
        ALLOCATE(eigenvalues(m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate eigenvalue array')
            DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK)
            RETURN
        END IF

        ALLOCATE(eigenvectors(n, m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate eigenvector array')
            DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, eigenvalues)
            RETURN
        END IF

        BLOCK
            REAL(rk), DIMENSION(m)    :: eval_real, eval_imag
            REAL(rk), DIMENSION(m, m) :: evec_right, evec_left
            REAL(rk), DIMENSION(:), ALLOCATABLE :: work_lapack
            INTEGER(ik) :: lwork_lapack, info_lapack
            REAL(rk)    :: work_query(1)

            CALL DGEEV('N', 'V', m, H_m, m, eval_real, eval_imag, &
                       evec_left, m, evec_right, m, work_query, -1, info_lapack)

            lwork_lapack = INT(work_query(1))
            ALLOCATE(work_lapack(lwork_lapack), STAT=ALLOC_STAT)
            IF (ALLOC_STAT /= 0) THEN
                ERROR_STATUS = ERR_ARNOLDI_ALLOC
                CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate LAPACK workspace')
                DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, eigenvalues, eigenvectors)
                RETURN
            END IF

            CALL DGEEV('N', 'V', m, H_m, m, eval_real, eval_imag, &
                       evec_left, m, evec_right, m, work_lapack, lwork_lapack, info_lapack)

            IF (info_lapack /= 0) THEN
                ERROR_STATUS = ERR_ARNOLDI_LAPACK
                CALL log_error(ERR_ARNOLDI_LAPACK, 'LAPACK DGEEV failed')
                DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, &
                           eigenvalues, eigenvectors, work_lapack)
                RETURN
            END IF

            DO i = 1, m
                eigenvalues(i) = CMPLX(eval_real(i), eval_imag(i), KIND=rk)
            END DO

            ! --- Step 5: Ritz vectors = V * (right eigenvectors of H_m) ---
            ! Handle real / complex-conjugate pairs as stored by DGEEV.
            i = 1
            DO WHILE (i <= m)
                IF (ABS(eval_imag(i)) < 1.0E-14_rk) THEN
                    CALL DGEMV('N', n, m, 1.0_rk, V, n, evec_right(:, i), 1, 0.0_rk, w, 1)
                    eigenvectors(:, i) = CMPLX(w, 0.0_rk, KIND=rk)
                    i = i + 1
                ELSE IF (eval_imag(i) > 0.0_rk) THEN
                    CALL DGEMV('N', n, m, 1.0_rk, V, n, evec_right(:, i),   1, 0.0_rk, w,      1)
                    CALL DGEMV('N', n, m, 1.0_rk, V, n, evec_right(:, i+1), 1, 0.0_rk, w_temp, 1)
                    eigenvectors(:, i) = CMPLX(w,  w_temp, KIND=rk)
                    IF (i+1 <= m) eigenvectors(:, i+1) = CMPLX(w, -w_temp, KIND=rk)
                    i = i + 2
                ELSE
                    i = i + 1
                END IF
            END DO

            DEALLOCATE(work_lapack)
        END BLOCK

        ! --- Step 6: Sort eigenvalues (and their Ritz vectors) ---
        ALLOCATE(sort_idx(m), imag_parts(m), temp_evec(n), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate sorting arrays')
            DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, eigenvalues, eigenvectors)
            RETURN
        END IF

        IF (PRESENT(sort_by)) THEN
            SELECT CASE (TRIM(sort_by))
            CASE ('imaginary')
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending imaginary part'
                DO i = 1, m
                    imag_parts(i) = AIMAG(eigenvalues(i))
                    sort_idx(i)   = i
                END DO
            CASE ('real')
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending real part'
                DO i = 1, m
                    imag_parts(i) = REAL(eigenvalues(i))
                    sort_idx(i)   = i
                END DO
            CASE DEFAULT
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending magnitude'
                DO i = 1, m
                    imag_parts(i) = ABS(eigenvalues(i))
                    sort_idx(i)   = i
                END DO
            END SELECT
        ELSE
            WRITE(*,'(A)') 'Sorting eigenvalues by: descending magnitude (default)'
            DO i = 1, m
                imag_parts(i) = ABS(eigenvalues(i))
                sort_idx(i)   = i
            END DO
        END IF

        ! Simple bubble sort (descending)
        DO i = 1, m-1
            DO j = i+1, m
                IF (imag_parts(sort_idx(j)) > imag_parts(sort_idx(i))) THEN
                    k           = sort_idx(i)
                    sort_idx(i) = sort_idx(j)
                    sort_idx(j) = k
                END IF
            END DO
        END DO

        DO i = 1, m
            eval_work(i)    = eigenvalues(sort_idx(i))
            evec_work(:, i) = eigenvectors(:, sort_idx(i))
        END DO

        eigenvalues  = eval_work
        eigenvectors = evec_work

        DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, &
                   sort_idx, imag_parts, temp_evec)

    END SUBROUTINE arnoldi_eigenvalues


    ! -------------------------------------------------------------------------
    ! apply_linearized_operator  --  STUB (TODO)
    !
    ! Intended behaviour:
    !     w = ( F(q0 + eps_0 * v) - F(q0) ) / eps_0
    ! where F is the nonlinear flow solver advanced by TTime (invoked via
    ! run_simulation in call_CFD). The resulting w approximates A*v, where
    ! A is the Jacobian of F linearized about the base state q0.
    !
    ! Suggested implementation outline:
    !   1. Retrieve / cache base state q0 (the current flow field).
    !   2. Form perturbed state q_pert = q0 + eps_0 * v_in, flatten-mapped
    !      back to (rho, U, V, W, p, T) and written to the OpenFOAM time dir.
    !   3. CALL run_simulation(COMMAND_RUN, ierr) to advance by TTime.
    !   4. Read the advanced perturbed state q_pert_adv.
    !   5. Ensure the advanced base state q0_adv is available (cache or
    !      re-advance q0 the same way).
    !   6. w_out = (q_pert_adv - q0_adv) / eps_0, in the same layout as v_in.
    !
    ! Until that's implemented, this stub sets w_out = 0 and returns
    ! ERR_ARNOLDI_NOT_IMPLEMENTED so the Arnoldi loop halts immediately.
    ! -------------------------------------------------------------------------
    SUBROUTINE apply_linearized_operator(v_in, w_out, frechet_order, eps_0, TTime, ierr)
        REAL(rk), DIMENSION(:), INTENT(IN)  :: v_in
        REAL(rk), DIMENSION(:), INTENT(OUT) :: w_out
        CHARACTER(len=*),       INTENT(IN)  :: frechet_order
        REAL(rk),               INTENT(IN)  :: eps_0
        REAL(rk),               INTENT(IN)  :: TTime
        INTEGER(ik),            INTENT(OUT) :: ierr

        ! NOTE: frechet_order, eps_0, TTime are INTENT(IN) placeholders for the
        !       real implementation and are intentionally unused in this stub.
        !       The compiler may emit "unused dummy argument" warnings; that is
        !       expected until the real matvec is implemented.

        w_out = 0.0_rk
        ierr  = ERR_ARNOLDI_NOT_IMPLEMENTED

        CALL log_error(ERR_ARNOLDI_NOT_IMPLEMENTED, &
            'apply_linearized_operator: Frechet matvec not yet implemented')

        ! Silence unused-argument warnings (no runtime effect).
        IF (.FALSE.) THEN
            w_out(1) = v_in(1) + eps_0 + TTime + REAL(LEN_TRIM(frechet_order), rk)
        END IF
    END SUBROUTINE apply_linearized_operator

END MODULE Arnoldi
