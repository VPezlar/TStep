MODULE read_flow
    USE accuracy
    USE variables
    USE setup
    USE OpenFOAM_IO

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

        ! Read configuration
        CALL configurationRead(error_status)
        IF (error_status /= 0) THEN
            WRITE(*,*) 'FATAL: Configuration read failed.'
            RETURN
        END IF

        WRITE(*,*) 'Attempting to read data from:', TRIM(file_grid)
        WRITE(*,*) 'Attempting to read data from:', TRIM(file_var)
        WRITE(*,*) 'Flowfield format:', TRIM(flow_format)

        IF (flow_format == 'OpenFOAM') THEN
            ! Read Flowfield Variables
            CALL read_OF_scalars(TRIM(file_var)//'p', N_HEADER_var, p_in, data_count, error_status)
            IF (error_status /= 0) THEN
                WRITE(*,*) 'FATAL: Failed to read pressure field.'
                RETURN
            END IF

            CALL read_OF_scalars(TRIM(file_var)//'rho', N_HEADER_var, rho_in, data_count, error_status)
            IF (error_status /= 0) THEN
                WRITE(*,*) 'FATAL: Failed to read density field.'
                RETURN
            END IF

            CALL read_OF_scalars(TRIM(file_var)//'T', N_HEADER_var, T_in, data_count, error_status)
            IF (error_status /= 0) THEN
                WRITE(*,*) 'FATAL: Failed to read temperature field.'
                RETURN
            END IF

            CALL read_OF_vectors(TRIM(file_var)//'U', N_HEADER_grid, U_in, V_in, W_in, data_count, error_status)
            IF (error_status /= 0) THEN
                WRITE(*,*) 'FATAL: Failed to read velocity field.'
                RETURN
            END IF

            ! Read Grid Coordinates (Cell Centers)
            CALL read_OF_vectors(TRIM(file_grid), N_HEADER_grid, Xgrid, Ygrid, Zgrid, data_count, error_status)
            IF (error_status /= 0) THEN
                WRITE(*,*) 'FATAL: Failed to read grid coordinates.'
                RETURN
            END IF

        END IF

        WRITE(*,*) 'SUCCESS: All data read successfully!'
        WRITE(*,*) 'Total data points:', data_count

    END SUBROUTINE read_flowfield

END MODULE read_flow
