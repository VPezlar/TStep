MODULE setup
    USE accuracy
    USE variables

    IMPLICIT NONE

CONTAINS

    SUBROUTINE configurationRead(ierr)
        INTEGER(ik), INTENT(out) :: ierr
        INTEGER(ik) :: unit_num, status_id
        CHARACTER(len=400) :: message

        NAMELIST / Setup / N_HEADER_grid, &
                           N_HEADER_var, &
                           file_grid, &
                           file_var, &
                           output_file, &
                           flow_format

        ierr = 0
        CALL get_unit(unit_num)
        OPEN(unit_num, file="../inputs/inputs.in", status="old", iostat=status_id, iomsg=message)

        IF (status_id /= 0) THEN
            WRITE(*,*) 'ERROR in configurationRead: Could not open file: ../inputs/inputs.in'
            WRITE(*,*) 'Message: ', TRIM(message)
            ierr = 1
            RETURN
        END IF

        READ(unit_num, nml=Setup, iostat=status_id, iomsg=message)

        IF (status_id /= 0) THEN
            WRITE(*,*) 'ERROR in configurationRead: Could not read namelist Setup'
            WRITE(*,*) 'Message: ', TRIM(message)
            ierr = 2
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
