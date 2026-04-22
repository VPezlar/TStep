! =============================================================================
! Arnoldi  --  Arnoldi iteration for the dominant Ritz eigenpairs of the
!              linearized flow operator A about the base state q0.
!
! Public:  arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime,
!                              rho0, p0, T0, U0, V0, W0,
!                              F_q0_vec, has_F_q0,
!                              eigenvalues, eigenvectors, ierr,
!                              skip_normalization, sort_by)
!   Builds a size-m Krylov basis V and Hessenberg H via modified Gram-Schmidt
!   with one reorthogonalization pass, computes eigenvalues/eigenvectors of
!   H_m with LAPACK DGEEV, lifts them to full-space Ritz vectors (V * y),
!   and sorts by 'magnitude' / 'real' / 'imaginary'.
!
! The Frechet matvec is done by apply_linearized_operator below:
!   w = sum_i weights(i) * F(q0 + alphas(i)*eps_0*v) / eps_0
! where (alphas, weights) come from frechet_stencil.Frechet_weights(order).
! For order = 1 the stencil is [0, 1] so the alpha=0 sample equals F(q0) --
! if the caller provides that cached (has_F_q0=.TRUE.), we reuse it and save
! one CFD solver call per Arnoldi iteration.
!
! See Mathias & Medeiros (2022), 'Optimal computational parameters for
! maximum accuracy and minimum cost of Arnoldi-based time-stepping methods'.
! =============================================================================
MODULE Arnoldi
    USE accuracy
    USE error_handling
    USE variables
    USE setup,          ONLY: stability_time_dir, clear_endpoint_folder, &
                              find_endpoint_folder
    USE read_flow,      ONLY: read_flowfield
    USE write_flow,     ONLY: write_flowfield
    USE call_CFD,       ONLY: run_simulation
    USE frechet_stencil,ONLY: Frechet_weights
    USE state_vector,   ONLY: NVARS, pack_state, unpack_state

    IMPLICIT NONE

    ! -------------------------------------------------------------------------
    ! Explicit interfaces for the LAPACK / BLAS routines we use.
    ! Declaring them here (rather than relying on implicit interfaces at the
    ! link step) lets the compiler catch argument-kind and rank mismatches at
    ! compile time. We match the standard reference LAPACK signatures using
    ! ik/rk (int32/real64), which is the native ABI of OpenBLAS/Netlib LAPACK
    ! on a typical 32-bit-integer build. If an ILP64 build is ever linked, the
    ! integer kind would need to change accordingly.
    ! -------------------------------------------------------------------------
    INTERFACE
        SUBROUTINE DGEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, &
                         VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
            IMPORT :: ik, rk
            CHARACTER(len=1), INTENT(IN)    :: JOBVL, JOBVR
            INTEGER(ik),      INTENT(IN)    :: N, LDA, LDVL, LDVR, LWORK
            INTEGER(ik),      INTENT(OUT)   :: INFO
            REAL(rk),         INTENT(INOUT) :: A(LDA, *)
            REAL(rk),         INTENT(OUT)   :: WR(*), WI(*)
            REAL(rk),         INTENT(OUT)   :: VL(LDVL, *), VR(LDVR, *)
            REAL(rk),         INTENT(OUT)   :: WORK(*)
        END SUBROUTINE DGEEV

        SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
            IMPORT :: ik, rk
            CHARACTER(len=1), INTENT(IN)    :: TRANS
            INTEGER(ik),      INTENT(IN)    :: M, N, LDA, INCX, INCY
            REAL(rk),         INTENT(IN)    :: ALPHA, BETA
            REAL(rk),         INTENT(IN)    :: A(LDA, *), X(*)
            REAL(rk),         INTENT(INOUT) :: Y(*)
        END SUBROUTINE DGEMV
    END INTERFACE

    PRIVATE
    PUBLIC :: arnoldi_eigenvalues

CONTAINS

    SUBROUTINE arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                                   rho0, p0, T0, U0, V0, W0, &
                                   F_q0_vec, has_F_q0, &
                                   eigenvalues, eigenvectors, ERROR_STATUS, &
                                   skip_normalization, sort_by)
        ! Computes Ritz eigenvalues and eigenvectors of the linearized flow
        ! operator using Arnoldi iteration.
        !
        ! The matvec w = A*v is supplied by the private routine
        ! apply_linearized_operator below, which approximates A via a
        ! Frechet-derivative finite difference using run_simulation.
        !
        ! Inputs:
        !   v_init      -- starting Krylov vector (length NVARS*Nmesh).
        !   m           -- requested Krylov subspace size (= krylov_size).
        !   frechet_order -- 1 or any positive even integer.
        !   eps_0, TTime -- Frechet magnitude and integration time tau.
        !   rho0..W0    -- the six per-cell arrays of the base state q0 held
        !                  in memory for the lifetime of the Arnoldi loop.
        !   F_q0_vec    -- F(q0) packed to length NVARS*Nmesh; only touched
        !                  when frechet_order=1 (alpha=0 sample reuse).
        !   has_F_q0    -- .TRUE. if F_q0_vec is populated. For order > 1
        !                  the alpha=0 node is never visited, so this flag
        !                  may safely be .FALSE. with a zero F_q0_vec.
        !
        ! If skip_normalization=.TRUE., the caller MUST ensure ||v_init|| = 1.

        ! Arguments
        REAL(rk), DIMENSION(:),              INTENT(IN)  :: v_init
        INTEGER(ik),                         INTENT(IN)  :: m
        INTEGER(ik),                         INTENT(IN)  :: frechet_order
        REAL(rk),                            INTENT(IN)  :: eps_0
        REAL(rk),                            INTENT(IN)  :: TTime
        REAL(rk), DIMENSION(:),              INTENT(IN)  :: rho0, p0, T0
        REAL(rk), DIMENSION(:),              INTENT(IN)  :: U0, V0, W0
        REAL(rk), DIMENSION(:),              INTENT(IN)  :: F_q0_vec
        LOGICAL,                             INTENT(IN)  :: has_F_q0
        COMPLEX(rk), DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: eigenvalues
        COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: eigenvectors
        INTEGER(ik),                         INTENT(OUT) :: ERROR_STATUS
        LOGICAL,          OPTIONAL,          INTENT(IN)  :: skip_normalization
        CHARACTER(len=*), OPTIONAL,          INTENT(IN)  :: sort_by

        ! Local Variables
        INTEGER(ik) :: n
        INTEGER(ik) :: m_eff   ! Effective Krylov size (may be < m on breakdown)
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
        n     = SIZE(v_init)
        m_eff = m   ! May shrink below if Arnoldi breakdown occurs

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
            WRITE(*,'(A,I0,A,I0,A)') '[Arnoldi] iteration ', j, ' of ', m, ': applying A'
            CALL apply_linearized_operator(V(:, j), w, frechet_order, eps_0, TTime, &
                                           rho0, p0, T0, U0, V0, W0, &
                                           F_q0_vec, has_F_q0, ERROR_STATUS)
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
                ! Truncate: only the leading j x j Hessenberg is reliable,
                ! so restrict DGEEV / sorting / output to m_eff = j.
                m_eff = j
                EXIT
            END IF
        END DO

        ! --- Step 3: Extract m×m upper block of Hessenberg matrix ---
        H_m = H(1:m, 1:m)

        ! --- Step 4: Eigenvalues and eigenvectors of H_m (LAPACK DGEEV) ---
        ALLOCATE(eigenvalues(m_eff), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate eigenvalue array')
            DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK)
            RETURN
        END IF

        ALLOCATE(eigenvectors(n, m_eff), STAT=ALLOC_STAT)
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

            CALL DGEEV('N', 'V', m_eff, H_m, m, eval_real, eval_imag, &
                       evec_left, m, evec_right, m, work_query, -1, info_lapack)

            lwork_lapack = INT(work_query(1))
            ALLOCATE(work_lapack(lwork_lapack), STAT=ALLOC_STAT)
            IF (ALLOC_STAT /= 0) THEN
                ERROR_STATUS = ERR_ARNOLDI_ALLOC
                CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate LAPACK workspace')
                DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, eigenvalues, eigenvectors)
                RETURN
            END IF

            CALL DGEEV('N', 'V', m_eff, H_m, m, eval_real, eval_imag, &
                       evec_left, m, evec_right, m, work_lapack, lwork_lapack, info_lapack)

            IF (info_lapack /= 0) THEN
                ERROR_STATUS = ERR_ARNOLDI_LAPACK
                CALL log_error(ERR_ARNOLDI_LAPACK, 'LAPACK DGEEV failed')
                DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, &
                           eigenvalues, eigenvectors, work_lapack)
                RETURN
            END IF

            !
            DO i = 1, m_eff
                eigenvalues(i) = CMPLX(eval_real(i), eval_imag(i), KIND=rk)
            END DO

            ! Keep eigenvalues as the raw DGEEV output mu_i -- these
            ! are the Ritz eigenvalues of the time-tau flow operator
            ! exp(tau*A). Do NOT apply LOG(mu)/tau here: LOG's principal
            ! branch gives artificial imaginary parts +- pi/tau for any
            ! negative-real mu (common for decaying, non-oscillatory
            ! modes), which polluted the results. Continuous-time
            ! lambda = log(mu)/tau is computed (with branch warnings)
            ! in write_eigen_files alongside mu itself.

            ! --- Step 5: Ritz vectors = V * (right eigenvectors of H_m) ---
            ! Handle real / complex-conjugate pairs as stored by DGEEV.
            i = 1
            DO WHILE (i <= m_eff)
                IF (ABS(eval_imag(i)) < 1.0E-14_rk) THEN
                    CALL DGEMV('N', n, m_eff, 1.0_rk, V, n, evec_right(:, i), 1, 0.0_rk, w, 1)
                    eigenvectors(:, i) = CMPLX(w, 0.0_rk, KIND=rk)
                    i = i + 1
                ELSE IF (eval_imag(i) > 0.0_rk) THEN
                    CALL DGEMV('N', n, m_eff, 1.0_rk, V, n, evec_right(:, i),   1, 0.0_rk, w,      1)
                    CALL DGEMV('N', n, m_eff, 1.0_rk, V, n, evec_right(:, i+1), 1, 0.0_rk, w_temp, 1)
                    eigenvectors(:, i) = CMPLX(w,  w_temp, KIND=rk)
                    IF (i+1 <= m_eff) eigenvectors(:, i+1) = CMPLX(w, -w_temp, KIND=rk)
                    i = i + 2
                ELSE
                    i = i + 1
                END IF
            END DO

            DEALLOCATE(work_lapack)
        END BLOCK

        ! --- Step 6: Sort eigenvalues (and their Ritz vectors) ---
        ALLOCATE(sort_idx(m_eff), imag_parts(m_eff), temp_evec(n), STAT=ALLOC_STAT)
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
                DO i = 1, m_eff
                    imag_parts(i) = AIMAG(eigenvalues(i))
                    sort_idx(i)   = i
                END DO
            CASE ('real')
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending real part'
                DO i = 1, m_eff
                    imag_parts(i) = REAL(eigenvalues(i))
                    sort_idx(i)   = i
                END DO
            CASE DEFAULT
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending magnitude'
                DO i = 1, m_eff
                    imag_parts(i) = ABS(eigenvalues(i))
                    sort_idx(i)   = i
                END DO
            END SELECT
        ELSE
            WRITE(*,'(A)') 'Sorting eigenvalues by: descending magnitude (default)'
            DO i = 1, m_eff
                imag_parts(i) = ABS(eigenvalues(i))
                sort_idx(i)   = i
            END DO
        END IF

        ! Simple bubble sort (descending)
        DO i = 1, m_eff-1
            DO j = i+1, m_eff
                IF (imag_parts(sort_idx(j)) > imag_parts(sort_idx(i))) THEN
                    k           = sort_idx(i)
                    sort_idx(i) = sort_idx(j)
                    sort_idx(j) = k
                END IF
            END DO
        END DO

        DO i = 1, m_eff
            eval_work(i)    = eigenvalues(sort_idx(i))
            evec_work(:, i) = eigenvectors(:, sort_idx(i))
        END DO

        eigenvalues(1:m_eff)     = eval_work(1:m_eff)
        eigenvectors(:, 1:m_eff) = evec_work(:, 1:m_eff)

        DEALLOCATE(V, H, H_m, w, w_temp, eval_work, evec_work, RWORK, &
                   sort_idx, imag_parts, temp_evec)

    END SUBROUTINE arnoldi_eigenvalues


    ! -------------------------------------------------------------------------
    ! apply_linearized_operator  --  Frechet matvec w = A * v.
    !
    ! Approximates the action of the Jacobian A = dF/dq|_{q0} on an arbitrary
    ! perturbation v via the finite-difference formula
    !     w = (1/eps_0) * sum_i weights(i) * F(q0 + alphas(i) * eps_0 * v)
    ! where (alphas, weights) is the Frechet stencil of the requested order.
    !
    ! Each non-zero stencil node requires one external CFD solver run:
    !   1. Build q_i (per-field) = q0 + alphas(i) * eps_0 * v.
    !   2. write_flowfield(q_i) -> <stability_dir>/1/ (overwrites the bodies;
    !      headers are preserved from the seed).
    !   3. run_simulation(COMMAND_RUN) advances <stab>/1/ to <stab>/<1+TTime>/.
    !   4. read_flowfield(<stab>/<1+TTime>/) -> F(q_i); pack into a flat
    !      vector and accumulate weights(i) * vec(F(q_i)) into w_out.
    !
    ! The alphas(i) == 0 node (hit only for frechet_order = 1) is handled
    ! specially: we reuse the cached F_q0_vec passed in by the caller
    ! (has_F_q0 == .TRUE.), skipping one CFD run per Arnoldi iteration.
    !
    ! Workspace arrays are allocated and freed per call (one matvec = one
    ! Arnoldi iteration); memory pressure is dominated by V, H which live
    ! in the driver routine.
    ! -------------------------------------------------------------------------
    SUBROUTINE apply_linearized_operator(v_in, w_out, frechet_order, eps_0, TTime, &
                                         rho0, p0, T0, U0, V0, W0, &
                                         F_q0_vec, has_F_q0, ierr)
        REAL(rk), DIMENSION(:), INTENT(IN)  :: v_in
        REAL(rk), DIMENSION(:), INTENT(OUT) :: w_out
        INTEGER(ik),            INTENT(IN)  :: frechet_order
        REAL(rk),               INTENT(IN)  :: eps_0
        REAL(rk),               INTENT(IN)  :: TTime  ! forwarded to solver via controlDict; unused here
        REAL(rk), DIMENSION(:), INTENT(IN)  :: rho0, p0, T0, U0, V0, W0
        REAL(rk), DIMENSION(:), INTENT(IN)  :: F_q0_vec
        LOGICAL,                INTENT(IN)  :: has_F_q0
        INTEGER(ik),            INTENT(OUT) :: ierr

        ! Stencil
        REAL(rk), DIMENSION(:), ALLOCATABLE :: alphas, weights

        ! Perturbation, unpacked per-field (length Nmesh each)
        REAL(rk), DIMENSION(:), ALLOCATABLE :: dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W

        ! Perturbed state q_i (length Nmesh each)
        REAL(rk), DIMENSION(:), ALLOCATABLE :: q_rho, q_p, q_T, q_U, q_V, q_W

        ! Solver output F(q_i) read back (length Nmesh each) + grid (scratch)
        REAL(rk), DIMENSION(:), ALLOCATABLE :: Fr, Fp, FT, FU, FV, FW
        REAL(rk), DIMENSION(:), ALLOCATABLE :: Xg, Yg, Zg
        REAL(rk), DIMENSION(:), ALLOCATABLE :: F_vec

        INTEGER(ik) :: Nmesh, data_count_read
        INTEGER(ik) :: i, alloc_stat
        CHARACTER(len=256) :: init_dir, end_dir

        ! TTime is an INTENT(IN) placeholder -- the actual integration time
        ! is controlled by <stability_dir>/system/controlDict. We suppress
        ! the unused-dummy warning with this harmless reference.
        IF (.FALSE.) THEN
            w_out(1) = w_out(1) + TTime
        END IF

        ierr   = 0
        Nmesh  = SIZE(rho0)
        init_dir = stability_time_dir(TSTEP_INITIAL_TIME)
        end_dir  = stability_time_dir(TSTEP_INITIAL_TIME + TTime)

        ! --- Sanity checks on in/out vector shapes ---
        IF (SIZE(v_in) /= NVARS * Nmesh .OR. SIZE(w_out) /= NVARS * Nmesh) THEN
            ierr = ERR_ARNOLDI_INVALID_DIM
            CALL log_error(ERR_ARNOLDI_INVALID_DIM, &
                'apply_linearized_operator: v_in/w_out length != NVARS*Nmesh')
            RETURN
        END IF

        ! --- Build Frechet stencil (nodes + weights) ---
        CALL Frechet_weights(frechet_order, alphas, weights, ierr)
        IF (ierr /= 0) RETURN

        ! --- Allocate workspaces (per-field) ---
        ALLOCATE(dv_rho(Nmesh), dv_p(Nmesh), dv_T(Nmesh), &
                 dv_U  (Nmesh), dv_V(Nmesh), dv_W(Nmesh), &
                 q_rho (Nmesh), q_p (Nmesh), q_T (Nmesh), &
                 q_U   (Nmesh), q_V (Nmesh), q_W (Nmesh), &
                 F_vec(NVARS*Nmesh), STAT=alloc_stat)
        IF (alloc_stat /= 0) THEN
            ierr = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, &
                'apply_linearized_operator: per-field workspace alloc failed')
            DEALLOCATE(alphas, weights)
            RETURN
        END IF

        ! --- Unpack v_in into its six per-field slices (read-only) ---
        CALL unpack_state(v_in, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, ierr)
        IF (ierr /= 0) THEN
            DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                       q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
            RETURN
        END IF

        ! --- Accumulate weights(i) * F(q0 + alphas(i)*eps_0*v) ---
        w_out = 0.0_rk
        DO i = 1, SIZE(alphas)

            IF (ABS(alphas(i)) < 1.0E-14_rk) THEN
                ! alpha == 0: reuse cached F(q0) if the caller provided it.
                ! Otherwise, fall through to a regular CFD run with q_i = q0.
                IF (has_F_q0) THEN
                    WRITE(*,'(A,ES12.4,A,ES12.4,A)') &
                        '  [matvec] alpha=', alphas(i), '  weight=', weights(i), &
                        '  (reusing cached F(q0))'
                    w_out = w_out + weights(i) * F_q0_vec
                    CYCLE
                END IF
            END IF

            ! Build perturbed state per-field: q = q0 + alpha*eps*dv.
            q_rho = rho0 + alphas(i) * eps_0 * dv_rho
            q_p   = p0   + alphas(i) * eps_0 * dv_p
            q_T   = T0   + alphas(i) * eps_0 * dv_T
            q_U   = U0   + alphas(i) * eps_0 * dv_U
            q_V   = V0   + alphas(i) * eps_0 * dv_V
            q_W   = W0   + alphas(i) * eps_0 * dv_W

            WRITE(*,'(A,ES12.4,A,ES12.4)') &
                '  [matvec] alpha=', alphas(i), '  weight=', weights(i)

            ! Write perturbed state into <stability_dir>/1/
            CALL write_flowfield(q_rho, q_p, q_T, q_U, q_V, q_W, Nmesh, ierr, &
                                 path_override=init_dir)
            IF (ierr /= 0) THEN
                DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                           q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
                RETURN
            END IF

            ! Clear stale endpoint folder, then advance the CFD solver by TTime.
            ! Critical: without this, OpenFOAM round-off (1.0 + 1000*1e-4 !=
            ! 1.1 exactly) collides with the existing 1.1/ folder and the
            ! solver auto-bumps timePrecision, writing its output elsewhere.
            ! TStep then reads stale data and Arnoldi breaks down immediately.
            CALL clear_endpoint_folder(ierr)
            CALL run_simulation(TRIM(COMMAND_RUN), ierr)
            IF (ierr /= 0) THEN
                DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                           q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
                RETURN
            END IF

            ! Robust post-solver endpoint lookup: ask setup for whichever
            ! numeric-named folder OpenFOAM actually wrote (1.1/,
            ! 1.099999999999989/, ...). Fully decouples us from
            ! timePrecision / IEEE-754 round-off in OpenFOAM's time loop.
            BLOCK
                CHARACTER(len=256) :: ep_mat
                LOGICAL :: ep_exists
                CALL find_endpoint_folder(ep_mat, ierr)
                IF (ierr /= 0) THEN
                    DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                               q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
                    RETURN
                END IF
                INQUIRE(FILE=TRIM(ep_mat)//'p', EXIST=ep_exists)
                IF (.NOT. ep_exists) THEN
                    ierr = ERR_ARNOLDI_INVALID_DIM
                    CALL log_error(ERR_ARNOLDI_INVALID_DIM, &
                        'matvec: endpoint folder '//TRIM(ep_mat)//&
                        ' has no p file.')
                    DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                               q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
                    RETURN
                END IF

                ! Read advanced state from the folder OpenFOAM just wrote.
                CALL read_flowfield(Fr, Fp, FT, FU, FV, FW, Xg, Yg, Zg, &
                                    data_count_read, ierr, path_override=ep_mat)
            END BLOCK
            IF (ierr /= 0) THEN
                DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                           q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
                IF (ALLOCATED(Fr)) DEALLOCATE(Fr, Fp, FT, FU, FV, FW, Xg, Yg, Zg)
                RETURN
            END IF

            IF (data_count_read /= Nmesh) THEN
                ierr = ERR_ARNOLDI_INVALID_DIM
                CALL log_error(ERR_ARNOLDI_INVALID_DIM, &
                    'matvec: post-solver data_count /= Nmesh')
                DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                           q_rho, q_p, q_T, q_U, q_V, q_W, F_vec, &
                           Fr, Fp, FT, FU, FV, FW, Xg, Yg, Zg)
                RETURN
            END IF

            ! Pack F(q_i) and accumulate.
            CALL pack_state(Fr, Fp, FT, FU, FV, FW, F_vec, ierr)
            IF (ierr /= 0) THEN
                DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                           q_rho, q_p, q_T, q_U, q_V, q_W, F_vec, &
                           Fr, Fp, FT, FU, FV, FW, Xg, Yg, Zg)
                RETURN
            END IF

            w_out = w_out + weights(i) * F_vec

            DEALLOCATE(Fr, Fp, FT, FU, FV, FW, Xg, Yg, Zg)
        END DO

        ! Divide by eps to complete the Frechet derivative.
        w_out = w_out / eps_0

        DEALLOCATE(alphas, weights, dv_rho, dv_p, dv_T, dv_U, dv_V, dv_W, &
                   q_rho, q_p, q_T, q_U, q_V, q_W, F_vec)
    END SUBROUTINE apply_linearized_operator

END MODULE Arnoldi
