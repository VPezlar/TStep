MODULE variables

    USE accuracy

    IMPLICIT NONE

    INTEGER(ik) :: N_HEADER_grid, N_HEADER_var
    CHARACTER(len=256) :: file_grid_in, file_var_in, file_var_out, output_file, flow_format, COMMAND_RUN
    REAL(rk) :: dist_mag

END MODULE variables
