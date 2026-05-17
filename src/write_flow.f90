! =============================================================================
! write_flow  --  writes flow fields BACK into the CFD solver's native format
!                 (currently OpenFOAM) so the solver can continue from them.
!
! Produces one file per field (p, rho, T, U) inside the target directory,
! which is an OpenFOAM time directory. No grid is written -- the mesh lives
! separately in the case's constant/polyMesh. Use this to hand perturbed or
! advanced flowfields back to rhoCentralFoam (or another solver).
!
! Path selection:
!   If path_override is PRESENT, writes scalars/vector there.
!   Otherwise writes to '<stability_dir>/1/' (composed via
!   setup.f90:stability_time_dir(TSTEP_INITIAL_TIME)), the default target
!   for a Frechet perturbation.
!
! IMPORTANT: the target directory must already contain OpenFOAM field files
! (p, rho, T, U) with valid headers. write_OF_scalars / write_OF_vectors in
! OpenFOAM_IO preserve the header of the existing target and only rewrite
! the numeric body. Use setup.f90:seed_stability_initial before the first
! call to populate <stability_dir>/1/ from <baseflow_field>.
!
! Contrast with write_output, which produces a single human-readable CSV for
! post-processing and analysis rather than for the solver.
! =============================================================================

MODULE write_flow
    USE accuracy
    USE variables
    USE setup
    USE OpenFOAM_IO
    USE SU2_IO
    USE error_handling

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: write_flowfield

CONTAINS

    SUBROUTINE write_flowfield(rho_out, p_out, T_out, U_out, V_out, W_out, &
                               data_count, error_status, path_override)
        REAL(rk), DIMENSION(:), INTENT(in) :: rho_out, p_out, T_out
        REAL(rk), DIMENSION(:), INTENT(in) :: U_out, V_out, W_out
        INTEGER(ik), INTENT(in) :: data_count
        INTEGER(ik), INTENT(out) :: error_status
        CHARACTER(len=*), OPTIONAL, INTENT(in) :: path_override

        CHARACTER(len=256) :: path_used

        error_status = 0

        ! Choose target: explicit override, or default based on flow_format.
        IF (PRESENT(path_override)) THEN
            path_used = path_override
        ELSE
            IF (TRIM(flow_format) == 'SU2') THEN
                path_used = su2_restart_in
            ELSE
                path_used = stability_time_dir(TSTEP_INITIAL_TIME)
            END IF
        END IF

        WRITE(*,*) 'Attempting to write data to:', TRIM(path_used)
        WRITE(*,*) 'Flowfield format:', TRIM(flow_format)

        IF (TRIM(flow_format) == 'OpenFOAM') THEN
            CALL write_OF_scalars(TRIM(path_used)//'p', N_HEADER_var, p_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_PRESSURE
                CALL log_error(ERR_FLOW_PRESSURE, 'File: '//TRIM(path_used)//'p')
                RETURN
            END IF

            CALL write_OF_scalars(TRIM(path_used)//'rho', N_HEADER_var, rho_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_DENSITY
                CALL log_error(ERR_FLOW_DENSITY, 'File: '//TRIM(path_used)//'rho')
                RETURN
            END IF

            CALL write_OF_scalars(TRIM(path_used)//'T', N_HEADER_var, T_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_TEMPERATURE
                CALL log_error(ERR_FLOW_TEMPERATURE, 'File: '//TRIM(path_used)//'T')
                RETURN
            END IF

            CALL write_OF_vectors(TRIM(path_used)//'U', N_HEADER_var, U_out, V_out, W_out, data_count, error_status)
            IF (error_status /= 0) THEN
                error_status = ERR_FLOW_VELOCITY
                CALL log_error(ERR_FLOW_VELOCITY, 'File: '//TRIM(path_used)//'U')
                RETURN
            END IF

        ELSE IF (TRIM(flow_format) == 'SU2') THEN
            CALL write_SU2_restart(TRIM(path_used), rho_out, p_out, T_out, &
                                  U_out, V_out, W_out, data_count, &
                                  error_status)
            IF (error_status /= 0) THEN
                CALL log_error(ERR_FLOW_PRESSURE, &
                    'SU2 write failed: '//TRIM(path_used))
                RETURN
            END IF

        ELSE
            error_status = ERR_FLOW_UNKNOWN_FORMAT
            CALL log_error(ERR_FLOW_UNKNOWN_FORMAT, &
                'flow_format = '//TRIM(flow_format)//' (supported: OpenFOAM, SU2)')
            RETURN
        END IF

        WRITE(*,*) 'SUCCESS: All data written successfully!'
        WRITE(*,*) 'Total data points:', data_count

    END SUBROUTINE write_flowfield

END MODULE write_flow      
