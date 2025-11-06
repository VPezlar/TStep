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
                              Xgrid, Ygrid, Zgrid, data_count, error_status)
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: rho_in, p_in, T_in
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: U_in, V_in, W_in
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: Xgrid, Ygrid, Zgrid
        INTEGER(ik), INTENT(out) :: data_count
        INTEGER(ik), INTENT(out) :: error_status

        error_status = 0

        WRITE(*,*) 'Attempting to read data from:', TRIM(file_grid_in)
        WRITE(*,*) 'Attempting to read data from:', TRIM(file_var_in)
        WRITE(*,*) 'Flowfield format:', TRIM(flow_format)

        IF (flow_format == 'OpenFOAM') THEN
            ! Read Flowfield Variables
            CALL read_OF_scalars(TRIM(file_var_in)//'p', N_HEADER_var, p_in, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_PRESSURE
                CALL log_error(ERR_FLOW_PRESSURE, 'File: '//TRIM(file_var_in)//'p')
                RETURN
            END IF

            CALL read_OF_scalars(TRIM(file_var_in)//'rho', N_HEADER_var, rho_in, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_DENSITY
                CALL log_error(ERR_FLOW_DENSITY, 'File: '//TRIM(file_var_in)//'rho')
                RETURN
            END IF

            CALL read_OF_scalars(TRIM(file_var_in)//'T', N_HEADER_var, T_in, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_TEMPERATURE
                CALL log_error(ERR_FLOW_TEMPERATURE, 'File: '//TRIM(file_var_in)//'T')
                RETURN
            END IF

            CALL read_OF_vectors(TRIM(file_var_in)//'U', N_HEADER_grid, U_in, V_in, W_in, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_VELOCITY
                CALL log_error(ERR_FLOW_VELOCITY, 'File: '//TRIM(file_var_in)//'U')
                RETURN
            END IF

            ! Read Grid Coordinates (Cell Centers)
            CALL read_OF_vectors(TRIM(file_grid_in), N_HEADER_grid, Xgrid, Ygrid, Zgrid, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_GRID
                CALL log_error(ERR_FLOW_GRID, 'File: '//TRIM(file_grid_in))
                RETURN
            END IF

        END IF

        WRITE(*,*) 'SUCCESS: All data read successfully!'
        WRITE(*,*) 'Total data points:', data_count

    END SUBROUTINE read_flowfield

END MODULE read_flow
