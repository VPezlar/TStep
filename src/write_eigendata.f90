! =============================================================================
! write_eigendata  --  writes Arnoldi Ritz eigenvalues and Ritz eigenvectors
!                      to plain-text .dat files for external analysis.
!
! Output convention (two files under <stability_dir>/output/):
!
!   eigenvalues.dat
!     Columns: idx  Re(mu)  Im(mu)  |mu|
!     mu is the DGEEV Ritz eigenvalue of the discrete time-tau operator
!     exp(tau*A). We deliberately DO NOT emit the continuous-time
!     eigenvalue lambda = log(mu)/tau anymore (see note below).
!
!   eigenvectors.dat (WIDE MATRIX FORMAT)
!     A flat, contiguous matrix designed for high-speed I/O.
!     The file contains a single header block, followed by Nmesh rows.
!     Each row contains the cell index and coordinates, followed by repeating
!     blocks of 12 columns for each computed mode:
!         idx  X  Y  Z  [Mode 1: Re(Rho) Im(Rho) Re(P) Im(P) ... Re(W) Im(W)] [Mode 2...]
!     This allows numpy.loadtxt() or readmatrix() to ingest the entire
!     data set into RAM in a single contiguous block without searching.
!
! Note: why log(mu)/tau was dropped
!     Earlier versions of this file emitted the continuous-time eigenvalue
!     lambda = log(mu)/tau alongside mu. The intrinsic Fortran LOG() is the
!     principal-branch complex logarithm, so whenever a Ritz value lands on
!     the negative real axis (Im(mu)=0, Re(mu)<0) LOG(mu) picks up a spurious
!     +i*pi, giving Im(lambda) = pi/tau. That is a numerical Nyquist-like
!     artefact of the branch cut, NOT a real frequency of the continuous
!     operator (the true imaginary part is ambiguous modulo 2*pi/tau in the
!     discrete-to-continuous map). We used to flag these points with an
!     OK/BC/ZERO "tag" column, but the cleaner fix is to just not perform
!     the transform here: mu is unambiguous, and any user who wants lambda
!     can compute log(mu)/tau in their own post-processing with the branch
!     convention appropriate to their problem (e.g. adding 2*pi*k/tau to
!     recover the correct Floquet harmonic). Dropping lambda also removes
!     the need for the tag column.
! =============================================================================
MODULE write_eigendata
    USE accuracy
    USE variables
    USE error_handling
    USE setup, ONLY: get_unit, stability_output_path, stability_time_dir, TSTEP_INITIAL_TIME

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: write_eigen_files

CONTAINS

    SUBROUTINE write_eigen_files(eigenvalues, eigenvectors, Xgrid, Ygrid, Zgrid, error_status)
        COMPLEX(rk), DIMENSION(:), INTENT(IN) :: eigenvalues
        COMPLEX(rk), DIMENSION(:,:), INTENT(IN) :: eigenvectors
        REAL(rk), DIMENSION(:), INTENT(IN) :: Xgrid, Ygrid, Zgrid
        INTEGER(ik), INTENT(OUT) :: error_status

        INTEGER(ik) :: unit_eval, unit_evec, io_stat
        INTEGER(ik) :: m_size, n_size, idx, jdx, n_write
        INTEGER(ik) :: Nmesh
        COMPLEX(rk) :: mu

        error_status = 0
        n_size = SIZE(eigenvectors, 1)
        m_size = SIZE(eigenvalues)

        ! Determine Nmesh locally from the 6 physical variables
        Nmesh = n_size / 6

        ! Sanity check grid size vs Arnoldi dimension
        IF (SIZE(Xgrid) /= Nmesh) THEN
            error_status = ERR_STATE_SHAPE_MISMATCH
            CALL log_error(ERR_STATE_SHAPE_MISMATCH, 'Grid size does not match Nmesh in eigendata')
            RETURN
        END IF

        ! =====================================================================
        ! eigenvalues.dat
        ! =====================================================================
        CALL get_unit(unit_eval)
        OPEN(UNIT=unit_eval, FILE=TRIM(stability_output_path('eigenvalues.dat')), &
             STATUS='REPLACE', ACTION='WRITE', IOSTAT=io_stat)

        IF (io_stat /= 0) THEN
            error_status = ERR_EIGENDATA_OPEN_EVAL
            CALL log_error(ERR_EIGENDATA_OPEN_EVAL, &
                'File: '//TRIM(stability_output_path('eigenvalues.dat')))
            RETURN
        END IF

        WRITE(unit_eval, '(A)') '# Ritz eigenvalues of exp(tau*A)  (DGEEV output mu)'
        WRITE(unit_eval, '(A,I0)')     '# Number of eigenvalues : ', m_size
        WRITE(unit_eval, '(A,A)')      '# Sorting criterion     : ', TRIM(eigenvalue_sort_by)
        WRITE(unit_eval, '(A,ES25.15)')'# Integration time tau  : ', TTime
        WRITE(unit_eval, '(A)') '#'
        WRITE(unit_eval, '(A)') '# The continuous-time eigenvalue lambda = log(mu)/tau is deliberately'
        WRITE(unit_eval, '(A)') '# NOT written: Fortran LOG() uses the principal branch, which creates'
        WRITE(unit_eval, '(A)') '# a spurious Im(lambda) = pi/tau whenever mu sits on the negative real'
        WRITE(unit_eval, '(A)') '# axis. Compute lambda in post-processing with the branch appropriate'
        WRITE(unit_eval, '(A)') '# to your physics (cf. Floquet harmonics k*2*pi/tau).'
        WRITE(unit_eval, '(A)') '#'
        WRITE(unit_eval, '(A)') '# Columns: idx   Re(mu)                    Im(mu)                    |mu|'

        DO idx = 1, m_size
            mu = eigenvalues(idx)
            WRITE(unit_eval, '(I6,3ES25.15)') idx, REAL(mu), AIMAG(mu), ABS(mu)
        END DO

        CLOSE(unit_eval)
        WRITE(*,*) 'SUCCESS: Wrote eigenvalues to ', TRIM(stability_output_path('eigenvalues.dat'))

        ! =====================================================================
        ! eigenvectors.dat  -- WIDE MATRIX FORMAT
        ! =====================================================================
        n_write = MIN(m_size, N_eig_write)

        CALL get_unit(unit_evec)
        OPEN(UNIT=unit_evec, FILE=TRIM(stability_output_path('eigenvectors.dat')), &
             STATUS='REPLACE', ACTION='WRITE', IOSTAT=io_stat)

        IF (io_stat /= 0) THEN
            error_status = ERR_EIGENDATA_OPEN_EVEC
            CALL log_error(ERR_EIGENDATA_OPEN_EVEC, &
                'File: '//TRIM(stability_output_path('eigenvectors.dat')))
            RETURN
        END IF

        ! Write unified header
        WRITE(unit_evec, '(A)')    '# Ritz eigenvectors of exp(tau*A) - WIDE MATRIX FORMAT'
        WRITE(unit_evec, '(A,I0)') '# Nmesh                       : ', Nmesh
        WRITE(unit_evec, '(A,I0)') '# System dimension n = 6*Nmesh: ', n_size
        WRITE(unit_evec, '(A,I0)') '# Number of Ritz modes computed : ', m_size
        WRITE(unit_evec, '(A,I0)') '# Number of Ritz modes on disk  : ', n_write
        WRITE(unit_evec, '(A,A)')  '# Sorting criterion             : ', TRIM(eigenvalue_sort_by)
        WRITE(unit_evec, '(A)') '#'
        WRITE(unit_evec, '(A)') '# Layout: One row per computational cell (Nmesh rows total).'
        WRITE(unit_evec, '(A)') '# Columns 1-4: idx, X, Y, Z'
        WRITE(unit_evec, '(A)') '# Subsequent columns are repeating blocks of 12 for each computed mode:'
        WRITE(unit_evec, '(A)') '#   [Re(Rho) Im(Rho) Re(P) Im(P) Re(T) Im(T) Re(U) Im(U) Re(V) Im(V) Re(W) Im(W)]'
        WRITE(unit_evec, '(A)') '#'

        ! Write the flat data matrix
        ! The 10000(1X,ES25.15) format descriptor safely allows up to 10000 floating point
        ! numbers on a single line, accommodating N_eig_write up to ~800.
        DO jdx = 1, Nmesh
            WRITE(unit_evec, '(I8, 3(1X,ES25.15), 10000(1X,ES25.15))') &
                jdx, Xgrid(jdx), Ygrid(jdx), Zgrid(jdx), &
                ( &
                  REAL (eigenvectors(          jdx, idx)), AIMAG(eigenvectors(          jdx, idx)), &
                  REAL (eigenvectors(  Nmesh + jdx, idx)), AIMAG(eigenvectors(  Nmesh + jdx, idx)), &
                  REAL (eigenvectors(2*Nmesh + jdx, idx)), AIMAG(eigenvectors(2*Nmesh + jdx, idx)), &
                  REAL (eigenvectors(3*Nmesh + jdx, idx)), AIMAG(eigenvectors(3*Nmesh + jdx, idx)), &
                  REAL (eigenvectors(4*Nmesh + jdx, idx)), AIMAG(eigenvectors(4*Nmesh + jdx, idx)), &
                  REAL (eigenvectors(5*Nmesh + jdx, idx)), AIMAG(eigenvectors(5*Nmesh + jdx, idx))  &
                , idx = 1, n_write )
        END DO

        CLOSE(unit_evec)
        WRITE(*,*) 'SUCCESS: Wrote eigenvectors to ', TRIM(stability_output_path('eigenvectors.dat'))
        WRITE(*,*)

    END SUBROUTINE write_eigen_files

END MODULE write_eigendata