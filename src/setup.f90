MODULE setup
    USE accuracy
    USE variables
    USE error_handling

    IMPLICIT NONE

CONTAINS

    SUBROUTINE configurationRead(ierr)
        INTEGER(ik), INTENT(out) :: ierr
        INTEGER(ik) :: unit_num, status_id
        CHARACTER(len=400) :: message

        NAMELIST / Setup / N_HEADER_grid, &
                           N_HEADER_var, &
                           file_grid_in, &
                           file_var_in, &
                           file_var_out, &
                           output_file, &
                           COMMAND_RUN, &
                           dist_mag, &
                           flow_format

        ierr = 0
        CALL get_unit(unit_num)
        OPEN(unit_num, file="../inputs/inputs.in", status="old", iostat=status_id, iomsg=message)

        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_FILE_OPEN
            CALL log_error(ERR_SETUP_FILE_OPEN, 'File: ../inputs/inputs.in - '//TRIM(message))
            RETURN
        END IF

        READ(unit_num, nml=Setup, iostat=status_id, iomsg=message)

        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_NAMELIST_READ
            CALL log_error(ERR_SETUP_NAMELIST_READ, TRIM(message))
            CLOSE(unit_num)
            RETURN
        END IF

        CLOSE(unit_num)

    END SUBROUTINE configurationRead

    ! Helper subroutine to safely get an unused file unit number.
    SUBROUTINE get_unit(u)
        INTEGER(ik) :: u
        INTEGER(ik) :: i
        LOGICAL :: is_opened

        DO i = 10, 99
            INQUIRE(unit=i, opened=is_opened)
            IF (.NOT. is_opened) THEN
                u = i
                RETURN
            END IF
        END DO

        u = 88
    END SUBROUTINE get_unit

END MODULE setup
