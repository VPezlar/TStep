! =============================================================================
! read_flow  --  high-level orchestration for loading a flowfield into
!                allocatable arrays in memory.
!
! `read_flowfield(rho,p,T,U,V,W, X,Y,Z, n, ierr [, path_override])` dispatches
! on flow_format (currently only 'OpenFOAM') and pulls five quantities from
! disk using the low-level readers in OpenFOAM_IO:
!   - scalars p, rho, T from <path>/{p,rho,T}
!   - vector  U (split into U,V,W) from <path>/U
!   - grid cell centers X,Y,Z from baseflow_grid (full path to C file)
!
! Path selection:
!   If path_override is PRESENT, reads scalars/vector from there.
!   Otherwise reads from `baseflow_field` (the archive, default).
!
! This lets the same routine serve two purposes:
!   - At startup: load U0 from baseflow_field (default behaviour)
!   - Inside Frechet loop: load solver output from the stability case's
!     end-time folder '<stability_dir>/<time_to_str(1+TTime)>/' (caller
!     passes path_override = setup.f90:stability_time_dir(1.0+TTime))
!
! Per-field data counts are verified against the first field read
! (ERR_FLOW_COUNT_MISMATCH) to catch inconsistent flowfield files early.
! Unknown flow_format values are rejected with ERR_FLOW_UNKNOWN_FORMAT.
!
! Contrast with read/write at the field-file level, which lives in OpenFOAM_IO.
! =============================================================================
MODULE read_flow
    USE accuracy
    USE variables
    USE setup
    USE OpenFOAM_IO
    USE error_handling

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_flowfield

CONTAINS

    SUBROUTINE read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, &
                              Xgrid, Ygrid, Zgrid, data_count, error_status, &
                              path_override)
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: rho_in, p_in, T_in
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: U_in, V_in, W_in
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: Xgrid, Ygrid, Zgrid
        INTEGER(ik), INTENT(out) :: data_count
        INTEGER(ik), INTENT(out) :: error_status
        CHARACTER(len=*), OPTIONAL, INTENT(in) :: path_override

        CHARACTER(len=256) :: path_used
        INTEGER(ik) :: count_p, count_rho, count_T, count_U, count_grid

        error_status = 0
        data_count   = 0

        ! Choose path: explicit override, or default baseflow archive
        IF (PRESENT(path_override)) THEN
            path_used = path_override
        ELSE
            path_used = baseflow_field
        END IF

        WRITE(*,*) 'Attempting to read data from:', TRIM(baseflow_grid)
        WRITE(*,*) 'Attempting to read data from:', TRIM(path_used)
        WRITE(*,*) 'Flowfield format:', TRIM(flow_format)

        IF (TRIM(flow_format) == 'OpenFOAM') THEN
            ! Pressure
            CALL read_OF_scalars(TRIM(path_used)//'p', N_HEADER_var, p_in, count_p, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_PRESSURE
                CALL log_error(ERR_FLOW_PRESSURE, 'File: '//TRIM(path_used)//'p')
                RETURN
            END IF
            data_count = count_p

            ! Density
            CALL read_OF_scalars(TRIM(path_used)//'rho', N_HEADER_var, rho_in, count_rho, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_DENSITY
                CALL log_error(ERR_FLOW_DENSITY, 'File: '//TRIM(path_used)//'rho')
                RETURN
            END IF
            IF (count_rho /= data_count) THEN
                error_status = ERR_FLOW_COUNT_MISMATCH
                CALL log_error(ERR_FLOW_COUNT_MISMATCH, &
                    'rho count '//TRIM(ADJUSTL(INT_TO_STR(count_rho)))// &
                    ' differs from p count '//TRIM(ADJUSTL(INT_TO_STR(data_count))))
                RETURN
            END IF

            ! Temperature
            CALL read_OF_scalars(TRIM(path_used)//'T', N_HEADER_var, T_in, count_T, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_TEMPERATURE
                CALL log_error(ERR_FLOW_TEMPERATURE, 'File: '//TRIM(path_used)//'T')
                RETURN
            END IF
            IF (count_T /= data_count) THEN
                error_status = ERR_FLOW_COUNT_MISMATCH
                CALL log_error(ERR_FLOW_COUNT_MISMATCH, &
                    'T count '//TRIM(ADJUSTL(INT_TO_STR(count_T)))// &
                    ' differs from p count '//TRIM(ADJUSTL(INT_TO_STR(data_count))))
                RETURN
            END IF

            ! Velocity (vector field; uses N_HEADER_var, not N_HEADER_grid)
            CALL read_OF_vectors(TRIM(path_used)//'U', N_HEADER_var, U_in, V_in, W_in, count_U, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_VELOCITY
                CALL log_error(ERR_FLOW_VELOCITY, 'File: '//TRIM(path_used)//'U')
                RETURN
            END IF
            IF (count_U /= data_count) THEN
                error_status = ERR_FLOW_COUNT_MISMATCH
                CALL log_error(ERR_FLOW_COUNT_MISMATCH, &
                    'U count '//TRIM(ADJUSTL(INT_TO_STR(count_U)))// &
                    ' differs from p count '//TRIM(ADJUSTL(INT_TO_STR(data_count))))
                RETURN
            END IF

            ! Grid cell centers (always from baseflow_grid, never overridden)
            CALL read_OF_vectors(TRIM(baseflow_grid), N_HEADER_grid, Xgrid, Ygrid, Zgrid, count_grid, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_GRID
                CALL log_error(ERR_FLOW_GRID, 'File: '//TRIM(baseflow_grid))
                RETURN
            END IF
            IF (count_grid /= data_count) THEN
                error_status = ERR_FLOW_COUNT_MISMATCH
                CALL log_error(ERR_FLOW_COUNT_MISMATCH, &
                    'grid count '//TRIM(ADJUSTL(INT_TO_STR(count_grid)))// &
                    ' differs from p count '//TRIM(ADJUSTL(INT_TO_STR(data_count))))
                RETURN
            END IF

        ELSE
            ! Unsupported format: fail loudly rather than return an empty flowfield.
            error_status = ERR_FLOW_UNKNOWN_FORMAT
            CALL log_error(ERR_FLOW_UNKNOWN_FORMAT, &
                'flow_format = '//TRIM(flow_format)//' (supported: OpenFOAM)')
            RETURN
        END IF

        WRITE(*,*) 'SUCCESS: All data read successfully!'
        WRITE(*,*) 'Total data points:', data_count

    END SUBROUTINE read_flowfield

END MODULE read_flow
