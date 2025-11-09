MODULE OpenFOAM_IO
    USE accuracy
    USE variables
    USE setup, ONLY: get_unit, INT_TO_STR
    USE error_handling

    IMPLICIT NONE
    ! Purpose: Provides subroutines to read custom datasets (Header, Count N, N lines of data).
    ! Error codes are now imported from error_handling module


CONTAINS

    ! The subroutine to read the column vector data (scalars) from the specified file.
    SUBROUTINE read_OF_scalars(filename, n_header_lines, data_vector, n_data_points, ierr)
        ! --- Subroutine Arguments ---
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER(ik), INTENT(in) :: n_header_lines
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: data_vector
        INTEGER(ik), INTENT(out) :: n_data_points
        INTEGER(ik), INTENT(out) :: ierr

        ! --- Local Variables ---
        INTEGER(ik) :: i
        INTEGER(ik) :: unit_num
        INTEGER(ik) :: count_read
        INTEGER(ik) :: iostat_val

        ! --- Initialization ---
        ierr = 0
        n_data_points = 0

        ! Safety check: Deallocate any previous data_vector array state before proceeding.
        IF (ALLOCATED(data_vector)) THEN
            DEALLOCATE(data_vector)
        END IF

        ! --- File Handling and Opening ---
        CALL get_unit(unit_num)

        ! Open the data file for reading.
        OPEN(unit=unit_num, file=TRIM(filename), status='old', &
            action='read', iostat=iostat_val)

        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_OPEN
            CALL log_error(ERR_SCALAR_OPEN, 'File: '//TRIM(filename))
            RETURN
        END IF

        ! 1. Skip Header Lines
        DO i = 1, n_header_lines
            READ(unit_num, *, iostat=iostat_val)
            IF (iostat_val /= 0) THEN
                ierr = ERR_SCALAR_HEADER
                CALL log_error(ERR_SCALAR_HEADER, 'File: '//TRIM(filename))
                CLOSE(unit_num)
                RETURN
            END IF
        END DO

        ! 2. Read the Number of Data Points (N)
        READ(unit_num, *, iostat=iostat_val) count_read
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_COUNT
            CALL log_error(ERR_SCALAR_COUNT, 'File: '//TRIM(filename))
            CLOSE(unit_num)
            RETURN
        END IF

        n_data_points = count_read

        ! Validation: Ensure the count is positive before allocation.
        IF (n_data_points <= 0) THEN
            ierr = ERR_SCALAR_INVALID_COUNT
            CALL log_error(ERR_SCALAR_INVALID_COUNT, 'File: '//TRIM(filename)//', Count: '//&
                          TRIM(ADJUSTL(INT_TO_STR(n_data_points))))
            CLOSE(unit_num)
            RETURN
        END IF

        ! 3. Skip the required single line after N (e.g., opening parenthesis)
        READ(unit_num, *, iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_SKIP
            CALL log_error(ERR_SCALAR_SKIP, 'File: '//TRIM(filename))
            CLOSE(unit_num)
            RETURN
        END IF

        ! 4. Allocate the Data Vector
        ALLOCATE(data_vector(n_data_points), stat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_ALLOC
            CALL log_error(ERR_SCALAR_ALLOC, 'File: '//TRIM(filename)//', Size: '//&
                          TRIM(ADJUSTL(INT_TO_STR(n_data_points))))
            CLOSE(unit_num)
            RETURN
        END IF

        ! 5. Read the Column Vector Data
        ! Reads all N values into the allocated array using free-format reading.
        READ(unit_num, *, iostat=iostat_val) data_vector(:)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_READ
            CALL log_error(ERR_SCALAR_READ, 'File: '//TRIM(filename))
        END IF

        ! Close the file unit, releasing it for future use.
        CLOSE(unit_num)

    END SUBROUTINE read_OF_scalars

    ! --- Subroutine 2: Reads X, Y, Z coordinate vectors ---
    SUBROUTINE read_OF_vectors(filename, n_header_lines, x_vector, y_vector, z_vector, n_data_points, ierr)
        ! --- Subroutine Arguments ---
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER(ik), INTENT(in) :: n_header_lines
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: x_vector
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: y_vector
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: z_vector
        INTEGER(ik), INTENT(out) :: n_data_points
        INTEGER(ik), INTENT(out) :: ierr

        ! --- Local Variables ---
        INTEGER(ik) :: i
        INTEGER(ik) :: unit_num
        INTEGER(ik) :: count_read
        INTEGER(ik) :: iostat_val
        CHARACTER(len=256) :: line_buffer
        INTEGER(ik) :: buf_len

        ! --- Initialization ---
        ierr = 0
        n_data_points = 0

        ! Safety check: Deallocate any previous arrays.
        IF (ALLOCATED(x_vector)) DEALLOCATE(x_vector)
        IF (ALLOCATED(y_vector)) DEALLOCATE(y_vector)
        IF (ALLOCATED(z_vector)) DEALLOCATE(z_vector)

        ! --- File Handling and Opening ---
        CALL get_unit(unit_num)

        ! Open the data file.
        OPEN(unit=unit_num, file=TRIM(filename), status='old', &
            action='read', iostat=iostat_val)

        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_OPEN
            CALL log_error(ERR_VECTOR_OPEN, 'File: '//TRIM(filename))
            RETURN
        END IF

        ! 1. Skip Header Lines
        DO i = 1, n_header_lines
            READ(unit_num, *, iostat=iostat_val)
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_HEADER
                CALL log_error(ERR_VECTOR_HEADER, 'File: '//TRIM(filename))
                CLOSE(unit_num)
                RETURN
            END IF
        END DO

        ! 2. Read the Number of Data Points (N)
        READ(unit_num, *, iostat=iostat_val) count_read
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_COUNT
            CALL log_error(ERR_VECTOR_COUNT, 'File: '//TRIM(filename))
            CLOSE(unit_num)
            RETURN
        END IF

        n_data_points = count_read

        ! 3. Skip the required single line after N
        READ(unit_num, *, iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_SKIP
            CALL log_error(ERR_VECTOR_SKIP, 'File: '//TRIM(filename))
            CLOSE(unit_num)
            RETURN
        END IF

        ! Validation and Allocation
        IF (n_data_points <= 0) THEN
            ierr = ERR_VECTOR_INVALID_COUNT
            CALL log_error(ERR_VECTOR_INVALID_COUNT, 'File: '//TRIM(filename)//', Count: '//&
                          TRIM(ADJUSTL(INT_TO_STR(n_data_points))))
            CLOSE(unit_num)
            RETURN
        END IF

        ! Allocate all three coordinate vectors simultaneously.
        ALLOCATE(x_vector(n_data_points), y_vector(n_data_points), z_vector(n_data_points), stat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_ALLOC
            CALL log_error(ERR_VECTOR_ALLOC, 'File: '//TRIM(filename)//', Size: '//&
                          TRIM(ADJUSTL(INT_TO_STR(n_data_points))))
            CLOSE(unit_num)
            RETURN
        END IF

        ! 4. Read the (X Y Z) Data (using the robust internal read loop)
        ! This loop is necessary to skip the mandatory external parentheses ( ) on each line.
        DO i = 1, n_data_points
            ! Read the entire line into a character buffer.
            READ(unit_num, '(A)', iostat=iostat_val) line_buffer
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_EOF
                CALL log_error(ERR_VECTOR_EOF, 'File: '//TRIM(filename)//', At index: '//&
                              TRIM(ADJUSTL(INT_TO_STR(i))))
                CLOSE(unit_num)
                RETURN
            END IF

            ! Find the length of the non-blank characters to strip off trailing blanks.
            buf_len = LEN_TRIM(line_buffer)

            ! Use an internal read to parse the numbers.
            ! The substring (2:buf_len-1) skips the opening '(' and closing ')'.
            READ(line_buffer(2:buf_len-1), *, iostat=iostat_val) x_vector(i), y_vector(i), z_vector(i)
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_FORMAT
                CALL log_error(ERR_VECTOR_FORMAT, 'File: '//TRIM(filename)//', Line: '//&
                              TRIM(ADJUSTL(INT_TO_STR(i))))
                CLOSE(unit_num)
                RETURN
            END IF
        END DO

        ! Close the file unit.
        CLOSE(unit_num)

    END SUBROUTINE read_OF_vectors

    ! Subroutine to write scalar data back to OpenFOAM file (preserving header and footer)
    SUBROUTINE write_OF_scalars(filename, n_header_lines, data_vector, n_data_points, ierr)
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER(ik), INTENT(in) :: n_header_lines
        REAL(rk), DIMENSION(:), INTENT(in) :: data_vector
        INTEGER(ik), INTENT(in) :: n_data_points
        INTEGER(ik), INTENT(out) :: ierr

        ! Local variables
        INTEGER(ik) :: i
        INTEGER(ik) :: unit_in, unit_out
        INTEGER(ik) :: iostat_val
        CHARACTER(len=256) :: line_buffer
        CHARACTER(len=256) :: temp_filename

        ierr = 0
        temp_filename = TRIM(filename)//'.tmp'

        ! Validation
        IF (n_data_points /= SIZE(data_vector)) THEN
            ierr = ERR_SCALAR_INVALID_COUNT
            CALL log_error(ERR_SCALAR_INVALID_COUNT, 'Mismatch between n_data_points and array size')
            RETURN
        END IF

        ! Get file unit for input and open it
        CALL get_unit(unit_in)

        ! Open original file for reading
        OPEN(unit=unit_in, file=TRIM(filename), status='old', action='read', iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_OPEN
            CALL log_error(ERR_SCALAR_OPEN, 'File: '//TRIM(filename))
            RETURN
        END IF

        ! Get file unit for output (after input is opened)
        CALL get_unit(unit_out)

        ! Open temporary file for writing
        OPEN(unit=unit_out, file=TRIM(temp_filename), status='replace', action='write', iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_OPEN
            CALL log_error(ERR_SCALAR_OPEN, 'File: '//TRIM(temp_filename))
            CLOSE(unit_in)
            RETURN
        END IF

        ! 1. Copy header lines
        DO i = 1, n_header_lines
            READ(unit_in, '(A)', iostat=iostat_val) line_buffer
            IF (iostat_val /= 0) THEN
                ierr = ERR_SCALAR_HEADER
                CALL log_error(ERR_SCALAR_HEADER, 'File: '//TRIM(filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
            WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
            IF (iostat_val /= 0) THEN
                ierr = ERR_SCALAR_HEADER
                CALL log_error(ERR_SCALAR_HEADER, 'Write failed: '//TRIM(temp_filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! 2. Copy the count line
        READ(unit_in, '(A)', iostat=iostat_val) line_buffer
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_COUNT
            CALL log_error(ERR_SCALAR_COUNT, 'File: '//TRIM(filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF
        WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_COUNT
            CALL log_error(ERR_SCALAR_COUNT, 'Write failed: '//TRIM(temp_filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF

        ! 3. Copy the opening parenthesis line
        READ(unit_in, '(A)', iostat=iostat_val) line_buffer
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_SKIP
            CALL log_error(ERR_SCALAR_SKIP, 'File: '//TRIM(filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF
        WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_SKIP
            CALL log_error(ERR_SCALAR_SKIP, 'Write failed: '//TRIM(temp_filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF

        ! 4. Skip old data in input file
        DO i = 1, n_data_points
            READ(unit_in, *, iostat=iostat_val)
            IF (iostat_val /= 0) THEN
                ierr = ERR_SCALAR_READ
                CALL log_error(ERR_SCALAR_READ, 'File: '//TRIM(filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! 5. Write new data values
        DO i = 1, n_data_points
            WRITE(unit_out, '(ES23.15E3)', iostat=iostat_val) data_vector(i)
            IF (iostat_val /= 0) THEN
                ierr = ERR_SCALAR_READ
                CALL log_error(ERR_SCALAR_READ, 'Write data failed: '//TRIM(temp_filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! 6. Copy remaining lines (closing parenthesis and any footer)
        DO
            READ(unit_in, '(A)', iostat=iostat_val) line_buffer
            IF (iostat_val /= 0) EXIT  ! End of file
            WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
            IF (iostat_val /= 0) THEN
                ierr = ERR_SCALAR_READ
                CALL log_error(ERR_SCALAR_READ, 'Write footer failed: '//TRIM(temp_filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! Close files
        CLOSE(unit_in)
        CLOSE(unit_out)

        ! Replace original file with temporary file using system command
        CALL EXECUTE_COMMAND_LINE('mv "'//TRIM(temp_filename)//'" "'//TRIM(filename)//'"', &
                                   EXITSTAT=iostat_val, CMDSTAT=i)
        IF (i /= 0 .OR. iostat_val /= 0) THEN
            ierr = ERR_SCALAR_OPEN
            CALL log_error(ERR_SCALAR_OPEN, 'Failed to rename temp file to: '//TRIM(filename))
            RETURN
        END IF

    END SUBROUTINE write_OF_scalars

    ! Subroutine to write vector data back to OpenFOAM file (preserving header and footer)
    SUBROUTINE write_OF_vectors(filename, n_header_lines, x_vector, y_vector, z_vector, n_data_points, ierr)
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER(ik), INTENT(in) :: n_header_lines
        REAL(rk), DIMENSION(:), INTENT(in) :: x_vector, y_vector, z_vector
        INTEGER(ik), INTENT(in) :: n_data_points
        INTEGER(ik), INTENT(out) :: ierr

        ! Local variables
        INTEGER(ik) :: i
        INTEGER(ik) :: unit_in, unit_out
        INTEGER(ik) :: iostat_val
        CHARACTER(len=256) :: line_buffer
        CHARACTER(len=256) :: temp_filename

        ierr = 0
        temp_filename = TRIM(filename)//'.tmp'

        ! Validation
        IF (n_data_points /= SIZE(x_vector) .OR. &
            n_data_points /= SIZE(y_vector) .OR. &
            n_data_points /= SIZE(z_vector)) THEN
            ierr = ERR_VECTOR_INVALID_COUNT
            CALL log_error(ERR_VECTOR_INVALID_COUNT, 'Mismatch between n_data_points and array sizes')
            RETURN
        END IF

        ! Get file unit for input and open it
        CALL get_unit(unit_in)

        ! Open original file for reading
        OPEN(unit=unit_in, file=TRIM(filename), status='old', action='read', iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_OPEN
            CALL log_error(ERR_VECTOR_OPEN, 'File: '//TRIM(filename))
            RETURN
        END IF

        ! Get file unit for output (after input is opened)
        CALL get_unit(unit_out)

        ! Open temporary file for writing
        OPEN(unit=unit_out, file=TRIM(temp_filename), status='replace', action='write', iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_OPEN
            CALL log_error(ERR_VECTOR_OPEN, 'File: '//TRIM(temp_filename))
            CLOSE(unit_in)
            RETURN
        END IF

        ! 1. Copy header lines
        DO i = 1, n_header_lines
            READ(unit_in, '(A)', iostat=iostat_val) line_buffer
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_HEADER
                CALL log_error(ERR_VECTOR_HEADER, 'File: '//TRIM(filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
            WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_HEADER
                CALL log_error(ERR_VECTOR_HEADER, 'Write failed: '//TRIM(temp_filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! 2. Copy the count line
        READ(unit_in, '(A)', iostat=iostat_val) line_buffer
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_COUNT
            CALL log_error(ERR_VECTOR_COUNT, 'File: '//TRIM(filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF
        WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_COUNT
            CALL log_error(ERR_VECTOR_COUNT, 'Write failed: '//TRIM(temp_filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF

        ! 3. Copy the opening parenthesis line
        READ(unit_in, '(A)', iostat=iostat_val) line_buffer
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_SKIP
            CALL log_error(ERR_VECTOR_SKIP, 'File: '//TRIM(filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF
        WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
        IF (iostat_val /= 0) THEN
            ierr = ERR_VECTOR_SKIP
            CALL log_error(ERR_VECTOR_SKIP, 'Write failed: '//TRIM(temp_filename))
            CLOSE(unit_in)
            CLOSE(unit_out)
            RETURN
        END IF

        ! 4. Skip old data in input file
        DO i = 1, n_data_points
            READ(unit_in, *, iostat=iostat_val)
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_READ
                CALL log_error(ERR_VECTOR_READ, 'File: '//TRIM(filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! 5. Write new vector data with parentheses
        DO i = 1, n_data_points
            WRITE(unit_out, '(A,ES23.15E3,A,ES23.15E3,A,ES23.15E3,A)', iostat=iostat_val) &
                '(', x_vector(i), ' ', y_vector(i), ' ', z_vector(i), ')'
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_READ
                CALL log_error(ERR_VECTOR_READ, 'Write data failed: '//TRIM(temp_filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! 6. Copy remaining lines (closing parenthesis and any footer)
        DO
            READ(unit_in, '(A)', iostat=iostat_val) line_buffer
            IF (iostat_val /= 0) EXIT  ! End of file
            WRITE(unit_out, '(A)', iostat=iostat_val) TRIM(line_buffer)
            IF (iostat_val /= 0) THEN
                ierr = ERR_VECTOR_READ
                CALL log_error(ERR_VECTOR_READ, 'Write footer failed: '//TRIM(temp_filename))
                CLOSE(unit_in)
                CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! Close files
        CLOSE(unit_in)
        CLOSE(unit_out)

        ! Replace original file with temporary file using system command
        CALL EXECUTE_COMMAND_LINE('mv "'//TRIM(temp_filename)//'" "'//TRIM(filename)//'"', &
                                   EXITSTAT=iostat_val, CMDSTAT=i)
        IF (i /= 0 .OR. iostat_val /= 0) THEN
            ierr = ERR_VECTOR_OPEN
            CALL log_error(ERR_VECTOR_OPEN, 'Failed to rename temp file to: '//TRIM(filename))
            RETURN
        END IF

    END SUBROUTINE write_OF_vectors

END MODULE OpenFOAM_IO
