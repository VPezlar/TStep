MODULE write_flow
    USE accuracy
    USE variables
    USE setup
    USE OpenFOAM_IO
    USE error_handling

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: write_flowfield

CONTAINS

    SUBROUTINE write_flowfield(rho_out, p_out, T_out, U_out, V_out, W_out, &
                                   data_count, error_status)
        REAL(rk), DIMENSION(:), INTENT(in) :: rho_out, p_out, T_out
        REAL(rk), DIMENSION(:), INTENT(in) :: U_out, V_out, W_out
        INTEGER(ik), INTENT(in) :: data_count
        INTEGER(ik), INTENT(out) :: error_status

        error_status = 0

        WRITE(*,*) 'Attempting to write data to:', TRIM(file_var)
        WRITE(*,*) 'Flowfield format:', TRIM(flow_format)

        IF (flow_format == 'OpenFOAM') THEN
            ! Write Flowfield Variables
            CALL write_OF_scalars(TRIM(file_var)//'p', N_HEADER_var, p_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_PRESSURE
                CALL log_error(ERR_FLOW_PRESSURE, 'File: '//TRIM(file_var)//'p')
                RETURN
            END IF

            CALL write_OF_scalars(TRIM(file_var)//'rho', N_HEADER_var, rho_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_DENSITY
                CALL log_error(ERR_FLOW_DENSITY, 'File: '//TRIM(file_var)//'rho')
                RETURN
            END IF

            CALL write_OF_scalars(TRIM(file_var)//'T', N_HEADER_var, T_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_TEMPERATURE
                CALL log_error(ERR_FLOW_TEMPERATURE, 'File: '//TRIM(file_var)//'T')
                RETURN
            END IF

            CALL write_OF_vectors(TRIM(file_var)//'U', N_HEADER_grid, U_out, V_out, W_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_VELOCITY
                CALL log_error(ERR_FLOW_VELOCITY, 'File: '//TRIM(file_var)//'U')
                RETURN
            END IF

        END IF

        WRITE(*,*) 'SUCCESS: All data written successfully!'
        WRITE(*,*) 'Total data points:', data_count

    END SUBROUTINE write_flowfield

END MODULE write_flow
