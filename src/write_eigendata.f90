! =============================================================================
! write_eigendata  --  writes Arnoldi Ritz eigenvalues and Ritz eigenvectors
!                      to plain-text .dat files for external analysis, and
!                      exports specific modes directly to OpenFOAM formats.
!
! Produces two files in <stability_dir>/:
!   eigenvalues.dat   --  one row per eigenvalue: index, Re, Im, |lambda|
!   eigenvectors.dat  --  base mesh block, then block-sorted physical variables
!                         (Rho, P, T, U, V, W) limited by N_eig_write.
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

    SUBROUTINE write_eigen_files(eigenvalues, eigenvectors, Xgrid, Ygrid, Zgrid, error_status)
        COMPLEX(rk), DIMENSION(:), INTENT(IN) :: eigenvalues
        COMPLEX(rk), DIMENSION(:,:), INTENT(IN) :: eigenvectors
        REAL(rk), DIMENSION(:), INTENT(IN) :: Xgrid, Ygrid, Zgrid
        INTEGER(ik), INTENT(OUT) :: error_status

        INTEGER(ik) :: unit_eval, unit_evec, io_stat
        INTEGER(ik) :: m_size, n_size, idx, jdx
        INTEGER(ik) :: Nmesh

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
        WRITE(unit_evec, '(A,I0)') '# Number of eigenvectors written to disk: ', MIN(m_size, N_eig_write)
        WRITE(unit_evec, '(A,A)') '# Sorting criterion: ', TRIM(eigenvalue_sort_by)
        WRITE(unit_evec, '(A)') '# Format: Block separated physical variables'
        WRITE(unit_evec, '(A)') '#         Component_Index  Real_Part  Imaginary_Part'
        WRITE(unit_evec, '(A)') '#'

        ! --- 1. Base Grid Block (Written Once) ---
        WRITE(unit_evec, '(A)') '# =========================================================='
        WRITE(unit_evec, '(A)') '# BASE MESH COORDINATES (X, Y, Z)'
        WRITE(unit_evec, '(A)') '# =========================================================='
        DO jdx = 1, Nmesh
            WRITE(unit_evec, '(I6,3ES25.15)') jdx, Xgrid(jdx), Ygrid(jdx), Zgrid(jdx)
        END DO
        WRITE(unit_evec, '(A)') '#'

        ! --- 2. Eigenvector Blocks ---
        ! Data: Write each eigenvector up to N_eig_write limit
        DO idx = 1, MIN(m_size, N_eig_write)
            WRITE(unit_evec, '(A,I0)') '# =========================================================='
            WRITE(unit_evec, '(A,I0)') '# EIGENVECTOR ', idx
            WRITE(unit_evec, '(A,I0)') '# =========================================================='

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


    ! =========================================================================
    ! export_mode_to_foam  -- Maps a flat Arnoldi eigenvector back into
    !                         OpenFOAM native format for ParaView.
    !
    ! Automatically clones the base directory header and splits the complex
    ! mode into Real, Imaginary, and Magnitude physical fields.
    ! =========================================================================
    SUBROUTINE export_mode_to_foam(idx, evec, Nmesh, error_status)
        USE state_vector, ONLY: unpack_state, NVARS
        USE write_flow, ONLY: write_flowfield

        IMPLICIT NONE
        INTEGER(ik), INTENT(IN) :: idx, Nmesh
        COMPLEX(rk), DIMENSION(:), INTENT(IN) :: evec
        INTEGER(ik), INTENT(OUT) :: error_status

        REAL(rk), ALLOCATABLE :: v_real(:), v_imag(:), v_mag(:)
        REAL(rk), ALLOCATABLE :: rho(:), p(:), T(:), U(:), V(:), W(:)
        CHARACTER(len=256) :: base_dir, target_dir
        CHARACTER(len=512) :: cmd
        INTEGER(ik) :: n, t_fake

        error_status = 0
        n = NVARS * Nmesh

        ! Allocate flat arrays and physical arrays
        ALLOCATE(v_real(n), v_imag(n), v_mag(n))
        ALLOCATE(rho(Nmesh), p(Nmesh), T(Nmesh), U(Nmesh), V(Nmesh), W(Nmesh))

        ! Deconstruct complex eigenvector
        v_real = REAL(evec)
        v_imag = AIMAG(evec)
        v_mag  = ABS(evec)

        ! The base directory holding the CFD headers (e.g. <stability_dir>/1/)
        base_dir = stability_time_dir(TSTEP_INITIAL_TIME)

        ! --- 1. Export REAL Part ---
        t_fake = 1000 + (idx * 10) + 1  ! e.g., 1011
        target_dir = TRIM(stability_time_dir(REAL(t_fake, rk)))
        cmd = 'mkdir -p ' // TRIM(target_dir) // ' && cp -f ' // TRIM(base_dir) // 'p ' // &
              TRIM(base_dir) // 'rho ' // TRIM(base_dir) // 'T ' // TRIM(base_dir) // 'U ' // TRIM(target_dir) // '/'
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE.)
        CALL unpack_state(v_real, rho, p, T, U, V, W, error_status)
        CALL write_flowfield(rho, p, T, U, V, W, Nmesh, error_status, path_override=target_dir)

        ! --- 2. Export IMAGINARY Part ---
        t_fake = 1000 + (idx * 10) + 2  ! e.g., 1012
        target_dir = TRIM(stability_time_dir(REAL(t_fake, rk)))
        cmd = 'mkdir -p ' // TRIM(target_dir) // ' && cp -f ' // TRIM(base_dir) // 'p ' // &
              TRIM(base_dir) // 'rho ' // TRIM(base_dir) // 'T ' // TRIM(base_dir) // 'U ' // TRIM(target_dir) // '/'
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE.)
        CALL unpack_state(v_imag, rho, p, T, U, V, W, error_status)
        CALL write_flowfield(rho, p, T, U, V, W, Nmesh, error_status, path_override=target_dir)

        ! --- 3. Export MAGNITUDE ---
        t_fake = 1000 + (idx * 10) + 3  ! e.g., 1013
        target_dir = TRIM(stability_time_dir(REAL(t_fake, rk)))
        cmd = 'mkdir -p ' // TRIM(target_dir) // ' && cp -f ' // TRIM(base_dir) // 'p ' // &
              TRIM(base_dir) // 'rho ' // TRIM(base_dir) // 'T ' // TRIM(base_dir) // 'U ' // TRIM(target_dir) // '/'
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE.)
        CALL unpack_state(v_mag, rho, p, T, U, V, W, error_status)
        CALL write_flowfield(rho, p, T, U, V, W, Nmesh, error_status, path_override=target_dir)

        DEALLOCATE(v_real, v_imag, v_mag, rho, p, T, U, V, W)
    END SUBROUTINE export_mode_to_foam

END MODULE write_eigendata