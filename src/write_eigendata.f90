! =============================================================================
! write_eigendata  --  writes Arnoldi Ritz eigenvalues and Ritz eigenvectors
!                      to plain-text .dat files for external analysis.
!
! Produces two files in ../output/:
!   eigenvalues.dat   --  one row per eigenvalue: index, Re, Im, |lambda|
!   eigenvectors.dat  --  for each eigenvector, n lines: component, Re, Im
! Both files are commented with a header recording m, n, and the sort order
! (`eigenvalue_sort_by` from inputs.in).
!
! Intended for post-processing with numpy/MATLAB etc. The module is only
! wired back into main once the Frechet matvec stub in Arnoldi.f90 is real.
! =============================================================================
MODULE write_eigendata
    USE accuracy
    USE variables
    USE error_handling
    USE setup, ONLY: get_unit
    
    IMPLICIT NONE
    
    PRIVATE
    PUBLIC :: write_eigen_files
    
CONTAINS

    SUBROUTINE write_eigen_files(eigenvalues, eigenvectors, error_status)
        ! Writes eigenvalues and eigenvectors to files in ../output/ directory
        ! Files: eigenvalues.dat and eigenvectors.dat
        ! Format: Sorted according to eigenvalue_sort_by from configuration
        !
        ! Arguments:
        !   eigenvalues   - Complex array of eigenvalues (size m)
        !   eigenvectors  - Complex array of eigenvectors (size n×m)
        !   error_status  - 0 for success, non-zero for error
        
        COMPLEX(rk), DIMENSION(:), INTENT(IN) :: eigenvalues
        COMPLEX(rk), DIMENSION(:,:), INTENT(IN) :: eigenvectors
        INTEGER(ik), INTENT(OUT) :: error_status
        
        INTEGER(ik) :: unit_eval, unit_evec, io_stat
        INTEGER(ik) :: n_size, m_size, idx, jdx
        
        error_status = 0
        n_size = SIZE(eigenvectors, 1)  ! System dimension
        m_size = SIZE(eigenvalues)      ! Number of eigenvalues (= Krylov size)
        
        ! --- Write eigenvalues.dat ---
        CALL get_unit(unit_eval)
        OPEN(UNIT=unit_eval, FILE='../output/eigenvalues.dat', STATUS='REPLACE', &
             ACTION='WRITE', IOSTAT=io_stat)
        
        IF (io_stat /= 0) THEN
            error_status = ERR_EIGENDATA_OPEN_EVAL
            CALL log_error(ERR_EIGENDATA_OPEN_EVAL, 'Check that ../output/ directory exists')
            RETURN
        END IF
        
        ! Header
        WRITE(unit_eval, '(A)') '# Ritz Eigenvalues from Arnoldi Iteration'
        WRITE(unit_eval, '(A,I0)') '# Number of eigenvalues: ', m_size
        WRITE(unit_eval, '(A,A)') '# Sorting criterion: ', TRIM(eigenvalue_sort_by)
        WRITE(unit_eval, '(A)') '# Format: Index  Real_Part  Imaginary_Part  Magnitude'
        WRITE(unit_eval, '(A)') '#'
        
        ! Data
        DO idx = 1, m_size
            WRITE(unit_eval, '(I6,3ES25.15)') idx, &
                REAL(eigenvalues(idx)), &
                AIMAG(eigenvalues(idx)), &
                ABS(eigenvalues(idx))
        END DO
        
        CLOSE(unit_eval)
        WRITE(*,*) 'SUCCESS: Wrote eigenvalues to ../output/eigenvalues.dat'
        
        ! --- Write eigenvectors.dat ---
        CALL get_unit(unit_evec)
        OPEN(UNIT=unit_evec, FILE='../output/eigenvectors.dat', STATUS='REPLACE', &
             ACTION='WRITE', IOSTAT=io_stat)
        
        IF (io_stat /= 0) THEN
            error_status = ERR_EIGENDATA_OPEN_EVEC
            CALL log_error(ERR_EIGENDATA_OPEN_EVEC, 'Check that ../output/ directory exists')
            RETURN
        END IF
        
        ! Header
        WRITE(unit_evec, '(A)') '# Ritz Eigenvectors from Arnoldi Iteration'
        WRITE(unit_evec, '(A,I0)') '# System dimension (n): ', n_size
        WRITE(unit_evec, '(A,I0)') '# Number of eigenvectors (m): ', m_size
        WRITE(unit_evec, '(A,A)') '# Sorting criterion: ', TRIM(eigenvalue_sort_by)
        WRITE(unit_evec, '(A)') '# Format: Each eigenvector spans n lines'
        WRITE(unit_evec, '(A)') '#         Component_Index  Real_Part  Imaginary_Part'
        WRITE(unit_evec, '(A)') '#'
        
        ! Data: Write each eigenvector
        DO idx = 1, m_size
            WRITE(unit_evec, '(A,I0)') '# Eigenvector ', idx
            DO jdx = 1, n_size
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), &
                    AIMAG(eigenvectors(jdx, idx))
            END DO
            WRITE(unit_evec, '(A)') '#'
        END DO
        
        CLOSE(unit_evec)
        WRITE(*,*) 'SUCCESS: Wrote eigenvectors to ../output/eigenvectors.dat'
        WRITE(*,*)
        
    END SUBROUTINE write_eigen_files

END MODULE write_eigendata
