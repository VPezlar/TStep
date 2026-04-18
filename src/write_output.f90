! =============================================================================
! write_output  --  writes a single human-readable CSV of the flowfield for
!                   post-processing, plotting, and debugging.
!
! Produces one flat file at `output_file` (e.g. ../output/flowfield.csv) with
! 9 columns: X,Y,Z,rho,p,T,U,V,W (one row per cell, scientific notation).
! Grid coordinates are included because a CSV has no implicit mesh. Format is
! solver-agnostic -- anything readable by pandas/MATLAB/ParaView etc.
!
! Contrast with write_flow, which writes back into the CFD solver's native
! format (OpenFOAM time directory) so the solver can continue stepping.
! =============================================================================

MODULE write_output
    USE accuracy
    USE variables
    USE setup
    USE error_handling

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: write_flowfield_data

CONTAINS

    SUBROUTINE write_flowfield_data(Xgrid, Ygrid, Zgrid, rho_in, p_in, T_in, &
                                    U_in, V_in, W_in, data_count, error_status)
        REAL(rk), DIMENSION(:), INTENT(in) :: Xgrid, Ygrid, Zgrid
        REAL(rk), DIMENSION(:), INTENT(in) :: rho_in, p_in, T_in
        REAL(rk), DIMENSION(:), INTENT(in) :: U_in, V_in, W_in
        INTEGER(ik), INTENT(in) :: data_count
        INTEGER(ik), INTENT(out) :: error_status

        INTEGER(ik) :: unit_num, iostat_val, i

        error_status = 0

        ! Get a free file unit
        CALL get_unit(unit_num)

        ! Open the output file
        OPEN(unit=unit_num, file=TRIM(output_file), status='replace', &
            action='write', iostat=iostat_val)

        IF (iostat_val /= 0) THEN
            error_status = ERR_OUTPUT_FILE_OPEN
            CALL log_error(ERR_OUTPUT_FILE_OPEN, 'File: '//TRIM(output_file))
            RETURN
        END IF

        ! Write header line with coordinates first
        WRITE(unit_num, '(A)') 'X,Y,Z,rho,p,T,U,V,W'

        ! Write data in column format (coordinates first, then fields)
        DO i = 1, data_count
            WRITE(unit_num, '(9(ES23.15E3,:,","))') Xgrid(i), Ygrid(i), Zgrid(i), &
                                                      rho_in(i), p_in(i), T_in(i), &
                                                      U_in(i), V_in(i), W_in(i)
        END DO

        CLOSE(unit_num)

        WRITE(*,*) 'SUCCESS: Flowfield data written to', TRIM(output_file)

    END SUBROUTINE write_flowfield_data

END MODULE write_output
