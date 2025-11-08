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

        ! General namelist (solver-agnostic)
        NAMELIST / General / flow_format, &
                             output_file, &
                             dist_mag, &
                             COMMAND_RUN

        ! Arnoldi-specific namelist
        NAMELIST / Arnoldi / krylov_size, &
                             frechet_order, &
                             eigenvalue_sort_by, &
                             eps_0, &
                             TTime

        ! OpenFOAM-specific namelist
        NAMELIST / OpenFOAM / N_HEADER_grid, &
                              N_HEADER_var, &
                              file_grid_in, &
                              file_var_in, &
                              file_var_out

        ierr = 0
        CALL get_unit(unit_num)
        OPEN(unit_num, file="../inputs/inputs.in", status="old", iostat=status_id, iomsg=message)

        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_FILE_OPEN
            CALL log_error(ERR_SETUP_FILE_OPEN, 'File: ../inputs/inputs.in - '//TRIM(message))
            RETURN
        END IF

        ! Read General namelist
        READ(unit_num, nml=General, iostat=status_id, iomsg=message)
        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_NAMELIST_READ
            CALL log_error(ERR_SETUP_NAMELIST_READ, 'General namelist - '//TRIM(message))
            CLOSE(unit_num)
            RETURN
        END IF

        ! Read Arnoldi namelist
        READ(unit_num, nml=Arnoldi, iostat=status_id, iomsg=message)
        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_NAMELIST_READ
            CALL log_error(ERR_SETUP_NAMELIST_READ, 'Arnoldi namelist - '//TRIM(message))
            CLOSE(unit_num)
            RETURN
        END IF

        ! Read solver-specific namelist based on flow_format
        IF (TRIM(flow_format) == 'OpenFOAM') THEN
            READ(unit_num, nml=OpenFOAM, iostat=status_id, iomsg=message)
            IF (status_id /= 0) THEN
                ierr = ERR_SETUP_NAMELIST_READ
                CALL log_error(ERR_SETUP_NAMELIST_READ, 'OpenFOAM namelist - '//TRIM(message))
                CLOSE(unit_num)
                RETURN
            END IF
        ELSE
            ! Future solvers can be added here
            ! ELSE IF (TRIM(flow_format) == 'SU2') THEN
            !     READ(unit_num, nml=SU2, iostat=status_id, iomsg=message)
            !     ...
            ierr = ERR_SETUP_INVALID_PARAM
            CALL log_error(ERR_SETUP_INVALID_PARAM, 'Unsupported flow_format: '//TRIM(flow_format))
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
