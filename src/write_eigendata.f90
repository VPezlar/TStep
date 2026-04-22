! =============================================================================
! write_eigendata  --  writes Arnoldi Ritz eigenvalues and Ritz eigenvectors
!                      to plain-text .dat files for external analysis, and
!                      exports specific modes directly to OpenFOAM formats.
!
! Produces two files in <stability_dir>/:
!   eigenvalues.dat   --  one row per eigenvalue: index, Re, Im, |lambda|
!   eigenvectors.dat  --  block sorted by physical variable (Rho, P, T, U, V, W)
!                         limited by N_eig_write from inputs.in
! =============================================================================
MODULE write_eigendata
    USE accuracy
    USE variables
    USE error_handling
    USE setup, ONLY: get_unit, stability_output_path, stability_time_dir, TSTEP_INITIAL_TIME

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: write_eigen_files, export_mode_to_foam

CONTAINS

    SUBROUTINE write_eigen_files(eigenvalues, eigenvectors, error_status)
        COMPLEX(rk), DIMENSION(:), INTENT(IN) :: eigenvalues
        COMPLEX(rk), DIMENSION(:,:), INTENT(IN) :: eigenvectors
        INTEGER(ik), INTENT(OUT) :: error_status

        INTEGER(ik) :: unit_eval, unit_evec, io_stat
        INTEGER(ik) :: m_size, n_size, idx, jdx
        INTEGER(ik) :: Nmesh

        error_status = 0
        n_size = SIZE(eigenvectors, 1)
        m_size = SIZE(eigenvalues)

        ! Determine Nmesh locally from the 6 physical variables
        Nmesh = n_size / 6

        ! --- Write eigenvalues.dat ---
        CALL get_unit(unit_eval)
        OPEN(UNIT=unit_eval, FILE=TRIM(stability_output_path('eigenvalues.dat')), &
             STATUS='REPLACE', ACTION='WRITE', IOSTAT=io_stat)

        IF (io_stat /= 0) THEN
            error_status = ERR_EIGENDATA_OPEN_EVAL
            CALL log_error(ERR_EIGENDATA_OPEN_EVAL, &
                'File: '//TRIM(stability_output_path('eigenvalues.dat')))
            RETURN
        END IF

        WRITE(unit_eval, '(A)') '# Ritz Eigenvalues from Arnoldi Iteration'
        WRITE(unit_eval, '(A,I0)') '# Number of eigenvalues: ', m_size
        WRITE(unit_eval, '(A,A)') '# Sorting criterion: ', TRIM(eigenvalue_sort_by)
        WRITE(unit_eval, '(A)') '# Format: Index  Real_Part  Imaginary_Part  Magnitude'
        WRITE(unit_eval, '(A)') '#'

        DO idx = 1, m_size
            WRITE(unit_eval, '(I6,3ES25.15)') idx, &
                REAL(eigenvalues(idx)), &
                AIMAG(eigenvalues(idx)), &
                ABS(eigenvalues(idx))
        END DO

        CLOSE(unit_eval)
        WRITE(*,*) 'SUCCESS: Wrote eigenvalues to ', TRIM(stability_output_path('eigenvalues.dat'))

        ! --- Write eigenvectors.dat ---
        CALL get_unit(unit_evec)
        OPEN(UNIT=unit_evec, FILE=TRIM(stability_output_path('eigenvectors.dat')), &
             STATUS='REPLACE', ACTION='WRITE', IOSTAT=io_stat)

        IF (io_stat /= 0) THEN
            error_status = ERR_EIGENDATA_OPEN_EVEC
            CALL log_error(ERR_EIGENDATA_OPEN_EVEC, &
                'File: '//TRIM(stability_output_path('eigenvectors.dat')))
            RETURN
        END IF

        WRITE(unit_evec, '(A)') '# Ritz Eigenvectors from Arnoldi Iteration'
        WRITE(unit_evec, '(A,I0)') '# System dimension (n): ', n_size
        WRITE(unit_evec, '(A,I0)') '# Number of eigenvectors computed (m): ', m_size
        WRITE(unit_evec, '(A,A)') '# Sorting criterion: ', TRIM(eigenvalue_sort_by)
        WRITE(unit_evec, '(A)') '# Format: Block separated physical variables'
        WRITE(unit_evec, '(A)') '#         Component_Index  Real_Part  Imaginary_Part'
        WRITE(unit_evec, '(A)') '#'

        ! Data: Write each eigenvector up to N_eig_write limit
        DO idx = 1, MIN(m_size, N_eig_write)
            WRITE(unit_evec, '(A,I0)') '# Eigenvector ', idx

            ! --- Rho ---
            WRITE(unit_evec, '(A)') '# Variable: Rho (Density)'
            DO jdx = 1, Nmesh
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), AIMAG(eigenvectors(jdx, idx))
            END DO

            ! --- Pressure ---
            WRITE(unit_evec, '(A)') '# Variable: P (Pressure)'
            DO jdx = Nmesh + 1, 2*Nmesh
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), AIMAG(eigenvectors(jdx, idx))
            END DO

            ! --- Temperature ---
            WRITE(unit_evec, '(A)') '# Variable: T (Temperature)'
            DO jdx = 2*Nmesh + 1, 3*Nmesh
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), AIMAG(eigenvectors(jdx, idx))
            END DO

            ! --- U Velocity ---
            WRITE(unit_evec, '(A)') '# Variable: U (X-Velocity)'
            DO jdx = 3*Nmesh + 1, 4*Nmesh
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), AIMAG(eigenvectors(jdx, idx))
            END DO

            ! --- V Velocity ---
            WRITE(unit_evec, '(A)') '# Variable: V (Y-Velocity)'
            DO jdx = 4*Nmesh + 1, 5*Nmesh
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), AIMAG(eigenvectors(jdx, idx))
            END DO

            ! --- W Velocity ---
            WRITE(unit_evec, '(A)') '# Variable: W (Z-Velocity)'
            DO jdx = 5*Nmesh + 1, 6*Nmesh
                WRITE(unit_evec, '(I6,2ES25.15)') jdx, &
                    REAL(eigenvectors(jdx, idx)), AIMAG(eigenvectors(jdx, idx))
            END DO

            WRITE(unit_evec, '(A)') '#'
        END DO

        CLOSE(unit_evec)
        WRITE(*,*) 'SUCCESS: Wrote eigenvectors to ', TRIM(stability_output_path('eigenvectors.dat'))
        WRITE(*,*)

    END SUBROUTINE write_eigen_files

END MODULE write_eigendata