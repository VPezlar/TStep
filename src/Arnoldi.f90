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
        ! Computes Ritz eigenvalues and eigenvectors using Arnoldi iteration
        !
        ! This subroutine performs Arnoldi iteration to build a Krylov subspace
        ! and then computes the Ritz values (approximate eigenvalues) and Ritz
        ! vectors (approximate eigenvectors) of the operator A.
        !
        ! The algorithm:
        !   1. Build Krylov basis V and Hessenberg matrix H via Arnoldi iteration
        !   2. Extract m×m upper block H_m from H
        !   3. Compute eigenvalues and eigenvectors of H_m using LAPACK
        !   4. Transform eigenvectors back to full space: Ritz vectors = V * y
        !   5. Sort by descending imaginary part
        !
        ! NOTE: Currently uses hardcoded test matrix A. This will be replaced
        !       with external CFD solver calls for the linearized Navier-Stokes operator.
        !
        ! IMPORTANT: If skip_normalization=.TRUE., the user MUST ensure that
        !            ||v_init|| = 1. Failure to do so will cause numerical errors
        !            in the Arnoldi iteration and incorrect eigenvalues.
        
        ! Arguments
        REAL(rk), DIMENSION(:), INTENT(IN)              :: v_init          ! Initial disturbance vector (size n)
        INTEGER(ik),            INTENT(IN)              :: m               ! Krylov subspace size
        CHARACTER(len=*),       INTENT(IN)              :: frechet_order   ! Order of Frechet derivative (unused for now)
        REAL(rk),               INTENT(IN)              :: eps_0           ! Epsilon parameter (unused for now)
        REAL(rk),               INTENT(IN)              :: TTime           ! Time parameter (unused for now)
        COMPLEX(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: eigenvalues    ! Ritz values (size m)
        COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: eigenvectors  ! Ritz vectors (size n×m)
        INTEGER(ik),            INTENT(OUT)             :: ERROR_STATUS    ! 0 for success, non-zero for error
        LOGICAL,                INTENT(IN), OPTIONAL   :: skip_normalization ! If .TRUE., skip v_init normalization (default: .FALSE.)
        CHARACTER(len=*),       INTENT(IN), OPTIONAL   :: sort_by         ! Sort criterion: 'magnitude', 'imaginary', 'real' (default: 'magnitude')
        
        ! Local Variables
        INTEGER(ik) :: n                              ! System dimension (inferred from v_init)
        REAL(rk), DIMENSION(:,:), ALLOCATABLE :: A    ! Test matrix (n×n) - will be replaced with CFD calls
        REAL(rk), DIMENSION(:,:), ALLOCATABLE :: V    ! Krylov basis vectors (n×(m+1))
        REAL(rk), DIMENSION(:,:), ALLOCATABLE :: H    ! Hessenberg matrix ((m+1)×m)
        REAL(rk), DIMENSION(:,:), ALLOCATABLE :: H_m  ! Upper m×m block of H for eigenvalue computation
        REAL(rk), DIMENSION(:), ALLOCATABLE :: w      ! Work vector
        REAL(rk) :: norm_w                            ! Norm of work vector
        REAL(rk) :: h_correction                      ! Correction for reorthogonalization
        INTEGER(ik) :: i, j, k                        ! Loop counters
        INTEGER(ik) :: ALLOC_STAT                     ! Allocation status
        
        ! Test matrix construction variables
        REAL(rk), DIMENSION(:,:), ALLOCATABLE :: eigvecs, eigvecs_inv
        REAL(rk), DIMENSION(:), ALLOCATABLE :: eigvals_diag
        REAL(rk), DIMENSION(:), ALLOCATABLE :: tau_qr, work_qr
        INTEGER(ik) :: info_decomp, lwork_qr
        
        ! LAPACK variables for eigenvalue computation
        COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: eval_work    ! Eigenvalues workspace (size m)
        COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: evec_work  ! Eigenvectors workspace (n×m)
        REAL(rk), DIMENSION(:), ALLOCATABLE :: RWORK           ! Real workspace for LAPACK
        
        ! Sorting variables
        INTEGER(ik), DIMENSION(:), ALLOCATABLE :: sort_idx
        REAL(rk), DIMENSION(:), ALLOCATABLE :: imag_parts
        COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: temp_evec
        
        ! Random number generator variables
        INTEGER(ik) :: clock_seed
        INTEGER(ik), DIMENSION(:), ALLOCATABLE :: seed_array
        INTEGER(ik) :: seed_size
        
        ! --- Initialization ---
        ERROR_STATUS = 0
        n = SIZE(v_init)
        
        ! Validate inputs
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
        ALLOCATE(A(n, n), V(n, m+1), H(m+1, m), H_m(m, m), w(n), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate Arnoldi arrays')
            RETURN
        END IF
        
        ALLOCATE(eigvecs(n, n), eigvecs_inv(n, n), eigvals_diag(n), tau_qr(n), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate matrix construction arrays')
            DEALLOCATE(A, V, H, H_m, w)
            RETURN
        END IF
        
        ALLOCATE(eval_work(m), evec_work(n, m), RWORK(2*m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate workspace arrays')
            DEALLOCATE(A, V, H, H_m, w)
            RETURN
        END IF
        
        ! --- Create test matrix A with known eigenvalues ---
        ! Strategy: A = Q @ diag(eigvals) @ Q^T
        ! where Q is orthogonal (from QR decomposition of random matrix)
        ! This ensures excellent numerical conditioning
        
        ! Generate linearly spaced eigenvalues: 1, 2, 3, ..., n
        DO i = 1, n
            eigvals_diag(i) = REAL(i, rk)
        END DO
        
        ! Generate random matrix and compute QR decomposition
        ! This gives us an orthogonal matrix Q (much better conditioned than random)
        CALL RANDOM_NUMBER(eigvecs)
        eigvecs = eigvecs - 0.5_rk  ! Center around zero
        
        ! Compute QR decomposition: eigvecs = Q @ R
        ! Query optimal work size
        ALLOCATE(work_qr(1))
        CALL DGEQRF(n, n, eigvecs, n, tau_qr, work_qr, -1, info_decomp)
        IF (info_decomp /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_LAPACK
            CALL log_error(ERR_ARNOLDI_LAPACK, 'Failed QR query')
            DEALLOCATE(A, V, H, H_m, w, eigvecs, eigvecs_inv, eigvals_diag, tau_qr, work_qr)
            RETURN
        END IF
        
        lwork_qr = INT(work_qr(1))
        DEALLOCATE(work_qr)
        ALLOCATE(work_qr(lwork_qr))
        
        ! Perform QR decomposition
        CALL DGEQRF(n, n, eigvecs, n, tau_qr, work_qr, lwork_qr, info_decomp)
        IF (info_decomp /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_LAPACK
            CALL log_error(ERR_ARNOLDI_LAPACK, 'Failed QR decomposition')
            DEALLOCATE(A, V, H, H_m, w, eigvecs, eigvecs_inv, eigvals_diag, tau_qr, work_qr)
            RETURN
        END IF
        
        ! Query optimal work size for DORGQR
        DEALLOCATE(work_qr)
        ALLOCATE(work_qr(1))
        CALL DORGQR(n, n, n, eigvecs, n, tau_qr, work_qr, -1, info_decomp)
        lwork_qr = INT(work_qr(1))
        DEALLOCATE(work_qr)
        ALLOCATE(work_qr(lwork_qr))
        
        ! Extract Q matrix
        CALL DORGQR(n, n, n, eigvecs, n, tau_qr, work_qr, lwork_qr, info_decomp)
        IF (info_decomp /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_LAPACK
            CALL log_error(ERR_ARNOLDI_LAPACK, 'Failed to generate Q matrix')
            DEALLOCATE(A, V, H, H_m, w, eigvecs, eigvecs_inv, eigvals_diag, tau_qr, work_qr)
            RETURN
        END IF
        
        ! Compute A = Q @ diag(eigvals) @ Q^T
        ! Step 1: Compute Q @ diag(eigvals)
        DO j = 1, n
            DO i = 1, n
                A(i, j) = eigvecs(i, j) * eigvals_diag(j)
            END DO
        END DO
        
        ! Step 2: A = (Q @ diag) @ Q^T
        A = MATMUL(A, TRANSPOSE(eigvecs))
        
        WRITE(*,'(A,I0,A,I0)') 'Arnoldi: Test matrix A (', n, 'x', n, ') with known eigenvalues'
        WRITE(*,'(A,F0.1,A,F0.1)') '         λ range: ', eigvals_diag(1), ' to ', eigvals_diag(n)
        
        DEALLOCATE(eigvecs, eigvecs_inv, eigvals_diag, tau_qr, work_qr)
        
        ! --- Step 1: Initialize and normalize first Krylov vector ---
        ! CRITICAL: For diagonal test matrix, use RANDOM vector to explore all eigenspaces
        ! This is the standard approach in classic Arnoldi (ARPACK, MATLAB eigs, etc.)
        ! The CFD initial vector may have structure that prevents exploring all modes
        
        ! Initialize random seed with system clock for different results each run
        CALL RANDOM_SEED(SIZE=seed_size)
        ALLOCATE(seed_array(seed_size))
        CALL SYSTEM_CLOCK(COUNT=clock_seed)
        seed_array = clock_seed + 37 * [(i-1, i=1, seed_size)]  ! Mix seed with different offsets
        CALL RANDOM_SEED(PUT=seed_array)
        DEALLOCATE(seed_array)
        
        ! Generate random vector in [-0.5, 0.5] and normalize
        CALL RANDOM_NUMBER(V(:, 1))
        V(:, 1) = V(:, 1) - 0.5_rk  ! Center around zero
        norm_w = NORM2(V(:, 1))
        V(:, 1) = V(:, 1) / norm_w
        
        WRITE(*,'(A)') 'Arnoldi: Using random initial vector for diagonal test (classic approach)'
        WRITE(*,'(A,ES12.5)') '         ||v_random|| = ', NORM2(V(:, 1))
        
        H = 0.0_rk
        
        ! --- Step 2: Arnoldi iteration ---
        DO j = 1, m
            ! Compute w = A * v_j (matrix-vector product)
            ! NOTE: This will be replaced with CFD solver call
            w = MATMUL(A, V(:, j))
            
            ! Classical Gram-Schmidt orthogonalization (first pass)
            DO i = 1, j
                H(i, j) = DOT_PRODUCT(V(:, i), w)
                w = w - H(i, j) * V(:, i)
            END DO
            
            ! Reorthogonalization pass (critical for numerical stability with large m)
            ! This is the standard approach used in ARPACK and other production codes
            ! Cost: 2x orthogonalization, but essential for m > 20-30
            DO i = 1, j
                h_correction = DOT_PRODUCT(V(:, i), w)
                H(i, j) = H(i, j) + h_correction  ! Accumulate total projection
                w = w - h_correction * V(:, i)     ! Remove remaining component
            END DO
            
            ! Compute norm and normalize
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
        
        ! --- Step 4: Compute eigenvalues and eigenvectors using LAPACK ---
        ! Use DGEEV for eigenvalues of real matrix
        ! Allocate arrays for eigenvalues and eigenvectors
        ALLOCATE(eigenvalues(m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate eigenvalue array')
            DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK)
            RETURN
        END IF
        
        ALLOCATE(eigenvectors(n, m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate eigenvector array')
            DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK, eigenvalues)
            RETURN
        END IF
        
        ! Call LAPACK DGEEV to compute eigenvalues and eigenvectors
        BLOCK
            REAL(rk), DIMENSION(m) :: eval_real, eval_imag
            REAL(rk), DIMENSION(m, m) :: evec_right, evec_left
            REAL(rk), DIMENSION(:), ALLOCATABLE :: work_lapack
            INTEGER(ik) :: lwork_lapack, info_lapack
            REAL(rk) :: work_query(1)
            
            ! Query optimal workspace size
            CALL DGEEV('N', 'V', m, H_m, m, eval_real, eval_imag, &
                       evec_left, m, evec_right, m, work_query, -1, info_lapack)
            
            lwork_lapack = INT(work_query(1))
            ALLOCATE(work_lapack(lwork_lapack), STAT=ALLOC_STAT)
            IF (ALLOC_STAT /= 0) THEN
                ERROR_STATUS = ERR_ARNOLDI_ALLOC
                CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate LAPACK workspace')
                DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK, eigenvalues, eigenvectors)
                RETURN
            END IF
            
            ! Compute eigenvalues and right eigenvectors
            CALL DGEEV('N', 'V', m, H_m, m, eval_real, eval_imag, &
                       evec_left, m, evec_right, m, work_lapack, lwork_lapack, info_lapack)
            
            IF (info_lapack /= 0) THEN
                ERROR_STATUS = ERR_ARNOLDI_LAPACK
                CALL log_error(ERR_ARNOLDI_LAPACK, 'LAPACK DGEEV failed')
                DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK, eigenvalues, eigenvectors, work_lapack)
                RETURN
            END IF
            
            ! Store eigenvalues as complex numbers
            DO i = 1, m
                eigenvalues(i) = CMPLX(eval_real(i), eval_imag(i), KIND=rk)
            END DO
            
            ! --- Step 5: Compute Ritz vectors (transform back to full space) ---
            ! Ritz vectors: eigenvectors = V * evec_right
            ! V is n×m, evec_right is m×m, result is n×m
            ! CRITICAL: Handle complex eigenvectors correctly
            ! When eigenvalues are complex conjugate pairs, DGEEV stores them specially:
            ! - If eval_imag(j) > 0: evec_right(:,j) + i*evec_right(:,j+1)
            ! - If eval_imag(j) < 0: evec_right(:,j) - i*evec_right(:,j+1)
            ! - If eval_imag(j) = 0: evec_right(:,j) is real
            
            i = 1
            DO WHILE (i <= m)
                IF (ABS(eval_imag(i)) < 1.0E-14_rk) THEN
                    ! Real eigenvalue: real eigenvector
                    eigenvectors(:, i) = CMPLX(MATMUL(V(:, 1:m), evec_right(:, i)), 0.0_rk, KIND=rk)
                    i = i + 1
                ELSE IF (eval_imag(i) > 0.0_rk) THEN
                    ! Complex conjugate pair: (λ, λ*) with eigenvectors (v, v*)
                    ! evec_right(:,i) is real part, evec_right(:,i+1) is imag part
                    eigenvectors(:, i) = CMPLX(MATMUL(V(:, 1:m), evec_right(:, i)), &
                                               MATMUL(V(:, 1:m), evec_right(:, i+1)), KIND=rk)
                    IF (i+1 <= m) THEN
                        eigenvectors(:, i+1) = CMPLX(MATMUL(V(:, 1:m), evec_right(:, i)), &
                                                     -MATMUL(V(:, 1:m), evec_right(:, i+1)), KIND=rk)
                    END IF
                    i = i + 2
                ELSE
                    ! This is the conjugate (already handled)
                    i = i + 1
                END IF
            END DO
            
            DEALLOCATE(work_lapack)
            
        END BLOCK
        
        ! --- Step 6: Sort eigenvalues ---
        ALLOCATE(sort_idx(m), imag_parts(m), temp_evec(n), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate sorting arrays')
            DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK, eigenvalues, eigenvectors)
            RETURN
        END IF
        
        ! Determine sort criterion
        IF (PRESENT(sort_by)) THEN
            IF (TRIM(sort_by) == 'imaginary') THEN
                ! Sort by descending imaginary part
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending imaginary part'
                DO i = 1, m
                    imag_parts(i) = AIMAG(eigenvalues(i))
                    sort_idx(i) = i
                END DO
            ELSE IF (TRIM(sort_by) == 'real') THEN
                ! Sort by descending real part
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending real part'
                DO i = 1, m
                    imag_parts(i) = REAL(eigenvalues(i))
                    sort_idx(i) = i
                END DO
            ELSE
                ! Default: sort by magnitude
                WRITE(*,'(A)') 'Sorting eigenvalues by: descending magnitude'
                DO i = 1, m
                    imag_parts(i) = ABS(eigenvalues(i))
                    sort_idx(i) = i
                END DO
            END IF
        ELSE
            ! Default: sort by magnitude
            WRITE(*,'(A)') 'Sorting eigenvalues by: descending magnitude (default)'
            DO i = 1, m
                imag_parts(i) = ABS(eigenvalues(i))
                sort_idx(i) = i
            END DO
        END IF
        
        ! Simple bubble sort (descending order)
        DO i = 1, m-1
            DO j = i+1, m
                IF (imag_parts(sort_idx(j)) > imag_parts(sort_idx(i))) THEN
                    k = sort_idx(i)
                    sort_idx(i) = sort_idx(j)
                    sort_idx(j) = k
                END IF
            END DO
        END DO
        
        ! Reorder eigenvalues and eigenvectors
        DO i = 1, m
            eval_work(i) = eigenvalues(sort_idx(i))
            evec_work(:, i) = eigenvectors(:, sort_idx(i))
        END DO
        
        eigenvalues = eval_work
        eigenvectors = evec_work
        
        ! --- Cleanup ---
        DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK, sort_idx, imag_parts, temp_evec)
        
    END SUBROUTINE arnoldi_eigenvalues

END MODULE Arnoldi
