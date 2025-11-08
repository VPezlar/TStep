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
        INTEGER(ik) :: i, j, k                        ! Loop counters
        INTEGER(ik) :: ALLOC_STAT                     ! Allocation status
        
        ! LAPACK variables for eigenvalue computation
        COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: eval_work    ! Eigenvalues workspace (size m)
        COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: evec_work  ! Eigenvectors workspace (n×m)
        REAL(rk), DIMENSION(:), ALLOCATABLE :: RWORK           ! Real workspace for LAPACK
        INTEGER(ik) :: LWORK, INFO                             ! LAPACK parameters (unused)
        
        ! Sorting variables
        INTEGER(ik), DIMENSION(:), ALLOCATABLE :: sort_idx
        REAL(rk), DIMENSION(:), ALLOCATABLE :: imag_parts
        COMPLEX(rk) :: temp_eval
        COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: temp_evec
        
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
        
        ALLOCATE(eval_work(m), evec_work(n, m), RWORK(2*m), STAT=ALLOC_STAT)
        IF (ALLOC_STAT /= 0) THEN
            ERROR_STATUS = ERR_ARNOLDI_ALLOC
            CALL log_error(ERR_ARNOLDI_ALLOC, 'Failed to allocate workspace arrays')
            DEALLOCATE(A, V, H, H_m, w)
            RETURN
        END IF
        
        ! --- Create test matrix A with known eigenvalues ---
        ! Using a tridiagonal matrix with specific structure
        ! Eigenvalues will be: λ_k = 2*cos(k*π/(n+1)) for k=1,...,n
        A = 0.0_rk
        DO i = 1, n
            A(i, i) = 0.0_rk
            IF (i < n) A(i, i+1) = 1.0_rk
            IF (i > 1) A(i, i-1) = 1.0_rk
        END DO
        
        WRITE(*,'(A,I0,A,I0)') 'Arnoldi: Using test matrix A (', n, 'x', n, ') - tridiagonal structure'
        
        ! --- Step 1: Initialize and normalize first Krylov vector ---
        V(:, 1) = v_init
        
        ! Check if normalization should be skipped
        IF (PRESENT(skip_normalization)) THEN
            IF (skip_normalization) THEN
                ! User has requested to skip normalization
                ! WARNING: User must ensure ||v_init|| = 1
                WRITE(*,'(A)') 'WARNING: Skipping initial vector normalization.'
                WRITE(*,'(A)') '         User must ensure ||v_init|| = 1 for correct results.'
                norm_w = NORM2(V(:, 1))
                WRITE(*,'(A,ES15.6)') '         Current ||v_init|| = ', norm_w
                IF (ABS(norm_w - 1.0_rk) > 1.0E-6_rk) THEN
                    WRITE(*,'(A)') '         WARNING: ||v_init|| significantly differs from 1!'
                    WRITE(*,'(A)') '                  Results may be inaccurate.'
                END IF
            ELSE
                ! Normalize as usual
                norm_w = NORM2(V(:, 1))
                IF (norm_w < 1.0E-12_rk) THEN
                    ERROR_STATUS = ERR_ARNOLDI_ZERO_V1
                    CALL log_error(ERR_ARNOLDI_ZERO_V1, 'Initial vector has zero or near-zero norm')
                    DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK)
                    RETURN
                END IF
                V(:, 1) = V(:, 1) / norm_w
            END IF
        ELSE
            ! Default behavior: normalize
            norm_w = NORM2(V(:, 1))
            IF (norm_w < 1.0E-12_rk) THEN
                ERROR_STATUS = ERR_ARNOLDI_ZERO_V1
                CALL log_error(ERR_ARNOLDI_ZERO_V1, 'Initial vector has zero or near-zero norm')
                DEALLOCATE(A, V, H, H_m, w, eval_work, evec_work, RWORK)
                RETURN
            END IF
            V(:, 1) = V(:, 1) / norm_w
        END IF
        
        H = 0.0_rk
        
        ! --- Step 2: Arnoldi iteration ---
        DO j = 1, m
            ! Compute w = A * v_j (matrix-vector product)
            ! NOTE: This will be replaced with CFD solver call
            w = MATMUL(A, V(:, j))
            
            ! Gram-Schmidt orthogonalization
            DO i = 1, j
                H(i, j) = DOT_PRODUCT(V(:, i), w)
                w = w - H(i, j) * V(:, i)
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
            eigenvectors = MATMUL(V(:, 1:m), evec_right)
            
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
