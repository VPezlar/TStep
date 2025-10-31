MODULE OpenFOAM_IO
    USE accuracy
    USE variables
    USE setup, ONLY: get_unit

    IMPLICIT NONE
    ! Purpose: Provides subroutines to read custom datasets (Header, Count N, N lines of data).

    ! Error codes for read_OF_scalars
    INTEGER(ik), PARAMETER :: ERR_SCALAR_OPEN = 1
    INTEGER(ik), PARAMETER :: ERR_SCALAR_HEADER = 2
    INTEGER(ik), PARAMETER :: ERR_SCALAR_COUNT = 3
    INTEGER(ik), PARAMETER :: ERR_SCALAR_INVALID_COUNT = 4
    INTEGER(ik), PARAMETER :: ERR_SCALAR_ALLOC = 5
    INTEGER(ik), PARAMETER :: ERR_SCALAR_SKIP = 6
    INTEGER(ik), PARAMETER :: ERR_SCALAR_READ = 7

    ! Error codes for read_OF_vectors
    INTEGER(ik), PARAMETER :: ERR_VECTOR_OPEN = 11
    INTEGER(ik), PARAMETER :: ERR_VECTOR_HEADER = 12
    INTEGER(ik), PARAMETER :: ERR_VECTOR_COUNT = 13
    INTEGER(ik), PARAMETER :: ERR_VECTOR_SKIP = 14
    INTEGER(ik), PARAMETER :: ERR_VECTOR_INVALID_COUNT = 15
    INTEGER(ik), PARAMETER :: ERR_VECTOR_ALLOC = 16
    INTEGER(ik), PARAMETER :: ERR_VECTOR_EOF = 17
    INTEGER(ik), PARAMETER :: ERR_VECTOR_FORMAT = 18


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
            WRITE(*,*) 'ERROR in read_OF_scalars: Could not open file:', TRIM(filename)
            ierr = ERR_SCALAR_OPEN
            RETURN
        END IF

        ! 1. Skip Header Lines
        DO i = 1, n_header_lines
            READ(unit_num, *, iostat=iostat_val)
            IF (iostat_val /= 0) THEN
                WRITE(*,*) 'ERROR in read_OF_scalars: Premature EOF while skipping header.'
                ierr = ERR_SCALAR_HEADER
                CLOSE(unit_num)
                RETURN
            END IF
        END DO

        ! 2. Read the Number of Data Points (N)
        READ(unit_num, *, iostat=iostat_val) count_read
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_scalars: Could not read data count (N).'
            ierr = ERR_SCALAR_COUNT
            CLOSE(unit_num)
            RETURN
        END IF

        n_data_points = count_read

        ! Validation: Ensure the count is positive before allocation.
        IF (n_data_points <= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_scalars: Invalid data count (N <= 0):', n_data_points
            ierr = ERR_SCALAR_INVALID_COUNT
            CLOSE(unit_num)
            RETURN
        END IF

        ! 3. Skip the required single line after N (e.g., opening parenthesis)
        READ(unit_num, *, iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_scalars: Premature EOF while skipping discard line.'
            ierr = ERR_SCALAR_SKIP
            CLOSE(unit_num)
            RETURN
        END IF

        ! 4. Allocate the Data Vector
        ALLOCATE(data_vector(n_data_points), stat=iostat_val)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_scalars: Could not allocate memory.'
            ierr = ERR_SCALAR_ALLOC
            CLOSE(unit_num)
            RETURN
        END IF

        ! 5. Read the Column Vector Data
        ! Reads all N values into the allocated array using free-format reading.
        READ(unit_num, *, iostat=iostat_val) data_vector(:)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'WARNING in read_OF_scalars: Could not read all data points.'
            ierr = ERR_SCALAR_READ
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
            WRITE(*,*) 'ERROR in read_OF_vectors: Could not open file:', TRIM(filename)
            ierr = ERR_VECTOR_OPEN
            RETURN
        END IF

        ! 1. Skip Header Lines
        DO i = 1, n_header_lines
            READ(unit_num, *, iostat=iostat_val)
            IF (iostat_val /= 0) THEN
                WRITE(*,*) 'ERROR in read_OF_vectors: Premature EOF while skipping header.'
                ierr = ERR_VECTOR_HEADER
                CLOSE(unit_num)
                RETURN
            END IF
        END DO

        ! 2. Read the Number of Data Points (N)
        READ(unit_num, *, iostat=iostat_val) count_read
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_vectors: Could not read data count (N).'
            ierr = ERR_VECTOR_COUNT
            CLOSE(unit_num)
            RETURN
        END IF

        n_data_points = count_read

        ! 3. Skip the required single line after N
        READ(unit_num, *, iostat=iostat_val)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_vectors: Premature EOF while skipping discard line.'
            ierr = ERR_VECTOR_SKIP
            CLOSE(unit_num)
            RETURN
        END IF

        ! Validation and Allocation
        IF (n_data_points <= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_vectors: Invalid data count (N <= 0):', n_data_points
            ierr = ERR_VECTOR_INVALID_COUNT
            CLOSE(unit_num)
            RETURN
        END IF

        ! Allocate all three coordinate vectors simultaneously.
        ALLOCATE(x_vector(n_data_points), y_vector(n_data_points), z_vector(n_data_points), stat=iostat_val)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) 'ERROR in read_OF_vectors: Could not allocate memory.'
            ierr = ERR_VECTOR_ALLOC
            CLOSE(unit_num)
            RETURN
        END IF

        ! 4. Read the (X Y Z) Data (using the robust internal read loop)
        ! This loop is necessary to skip the mandatory external parentheses ( ) on each line.
        DO i = 1, n_data_points
            ! Read the entire line into a character buffer.
            READ(unit_num, '(A)', iostat=iostat_val) line_buffer
            IF (iostat_val /= 0) THEN
                WRITE(*,*) 'WARNING in read_OF_vectors: EOF encountered prematurely at index', i
                ierr = ERR_VECTOR_EOF
                CLOSE(unit_num)
                RETURN
            END IF

            ! Find the length of the non-blank characters to strip off trailing blanks.
            buf_len = LEN_TRIM(line_buffer)

            ! Use an internal read to parse the numbers.
            ! The substring (2:buf_len-1) skips the opening '(' and closing ')'.
            READ(line_buffer(2:buf_len-1), *, iostat=iostat_val) x_vector(i), y_vector(i), z_vector(i)
            IF (iostat_val /= 0) THEN
                WRITE(*,*) 'ERROR in read_OF_vectors: Format error on line:', i
                ierr = ERR_VECTOR_FORMAT
                CLOSE(unit_num)
                RETURN
            END IF
        END DO

        ! Close the file unit.
        CLOSE(unit_num)

    END SUBROUTINE read_OF_vectors

END MODULE OpenFOAM_IO
