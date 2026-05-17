! =============================================================================
! SU2_IO  --  low-level readers and writers for SU2 ASCII restart files.
!
! SU2 restart format (OUTPUT_FILES=(RESTART_ASCII)):
!   Line 1:  comma-separated double-quoted column names (header)
!   Lines 2..N+1:  comma-separated numeric data rows
!   Trailing lines: optional metadata (non-numeric first token)
!
! Column indices are discovered at runtime by parsing the header line;
! no hard-coded positions. Required columns:
!   x, y, z, Density, Momentum_x, Momentum_y, Momentum_z, Energy
!
! All conservative <-> primitive conversion happens here:
!   Conservative (SU2): rho, rho*u, rho*v, rho*w, rho*E
!   Primitive  (TStep): rho, p, T, U, V, W
!
! Two public routines:
!   read_SU2_solution  -- read restart, convert to primitive, extract grid
!   write_SU2_restart  -- convert primitive to conservative, write restart
!                         preserving header and non-flow columns
!
! Styled after OpenFOAM_IO.f90 (same USE list, error_handling idiom,
! get_unit, c_rename pattern).
! =============================================================================
MODULE SU2_IO
    USE accuracy
    USE variables
    USE setup, ONLY: get_unit, INT_TO_STR
    USE error_handling
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR

    IMPLICIT NONE

    ! libc rename(3) for atomic temp->final file replacement
    INTERFACE
        FUNCTION c_rename(oldpath, newpath) BIND(C, NAME="rename") RESULT(r)
            IMPORT :: C_INT, C_CHAR
            CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: oldpath, newpath
            INTEGER(C_INT) :: r
        END FUNCTION c_rename
    END INTERFACE

    ! Maximum number of columns we expect in a SU2 restart
    INTEGER(ik), PARAMETER :: MAX_SU2_COLS = 200

    PRIVATE
    PUBLIC :: read_SU2_solution, write_SU2_restart

CONTAINS

    ! -------------------------------------------------------------------------
    ! parse_header  --  split a comma-separated header of quoted column names
    !                   into an array of unquoted name strings.
    !
    ! Input:  header_line  e.g. '"PointID","x","y","z","Density",...'
    ! Output: names(1:ncols) with quotes stripped, ncols = column count
    ! -------------------------------------------------------------------------
    SUBROUTINE parse_header(header_line, names, ncols)
        CHARACTER(len=*), INTENT(IN)  :: header_line
        CHARACTER(len=64), INTENT(OUT) :: names(MAX_SU2_COLS)
        INTEGER(ik), INTENT(OUT) :: ncols

        INTEGER :: pos, start, hlen, j
        CHARACTER(len=64) :: token

        ncols = 0
        hlen  = LEN_TRIM(header_line)
        pos   = 1

        DO WHILE (pos <= hlen)
            ! Skip leading whitespace and commas
            DO WHILE (pos <= hlen .AND. &
                      (header_line(pos:pos) == ',' .OR. &
                       header_line(pos:pos) == ' '))
                pos = pos + 1
            END DO
            IF (pos > hlen) EXIT

            ! Find the token: from pos to next comma or end
            start = pos
            DO WHILE (pos <= hlen .AND. header_line(pos:pos) /= ',')
                pos = pos + 1
            END DO

            token = ADJUSTL(header_line(start:pos-1))

            ! Strip surrounding double quotes
            j = LEN_TRIM(token)
            IF (j >= 2) THEN
                IF (token(1:1) == '"' .AND. token(j:j) == '"') THEN
                    token = token(2:j-1)
                END IF
            END IF

            ncols = ncols + 1
            IF (ncols <= MAX_SU2_COLS) THEN
                names(ncols) = TRIM(token)
            END IF
        END DO
    END SUBROUTINE parse_header


    ! -------------------------------------------------------------------------
    ! find_column  --  locate a column name in the parsed header array.
    !                  Returns the 1-based index or 0 if not found.
    ! -------------------------------------------------------------------------
    FUNCTION find_column(names, ncols, target_name) RESULT(idx)
        CHARACTER(len=64), INTENT(IN) :: names(MAX_SU2_COLS)
        INTEGER(ik), INTENT(IN)  :: ncols
        CHARACTER(len=*), INTENT(IN) :: target_name
        INTEGER(ik) :: idx
        INTEGER(ik) :: i

        idx = 0
        DO i = 1, ncols
            IF (TRIM(names(i)) == TRIM(target_name)) THEN
                idx = i
                RETURN
            END IF
        END DO
    END FUNCTION find_column


    ! -------------------------------------------------------------------------
    ! is_data_row  --  returns .TRUE. if the first comma-separated token of
    !                  `line` parses as a number (integer or real).
    ! -------------------------------------------------------------------------
    FUNCTION is_data_row(line) RESULT(is_data)
        CHARACTER(len=*), INTENT(IN) :: line
        LOGICAL :: is_data
        CHARACTER(len=64) :: token
        INTEGER :: comma_pos, ios
        REAL(rk) :: dummy

        is_data = .FALSE.
        IF (LEN_TRIM(line) == 0) RETURN

        ! Extract first token (before first comma)
        comma_pos = INDEX(line, ',')
        IF (comma_pos > 1) THEN
            token = ADJUSTL(line(1:comma_pos-1))
        ELSE IF (comma_pos == 0) THEN
            token = ADJUSTL(TRIM(line))
        ELSE
            RETURN
        END IF

        ! Try to parse as a number
        READ(token, *, IOSTAT=ios) dummy
        is_data = (ios == 0)
    END FUNCTION is_data_row


    ! -------------------------------------------------------------------------
    ! parse_csv_row  --  extract the value at column `col_idx` from a
    !                    comma-separated line. Returns it as a real.
    ! -------------------------------------------------------------------------
    SUBROUTINE parse_csv_value(line, col_idx, val, ierr)
        CHARACTER(len=*), INTENT(IN) :: line
        INTEGER(ik), INTENT(IN)  :: col_idx
        REAL(rk), INTENT(OUT) :: val
        INTEGER(ik), INTENT(OUT) :: ierr

        INTEGER :: pos, start, llen, cur_col, ios

        ierr    = 0
        val     = 0.0_rk
        llen    = LEN_TRIM(line)
        pos     = 1
        cur_col = 0

        DO WHILE (pos <= llen)
            ! Skip leading spaces
            DO WHILE (pos <= llen .AND. line(pos:pos) == ' ')
                pos = pos + 1
            END DO

            start = pos
            ! Advance to next comma or end
            DO WHILE (pos <= llen .AND. line(pos:pos) /= ',')
                pos = pos + 1
            END DO

            cur_col = cur_col + 1
            IF (cur_col == col_idx) THEN
                READ(line(start:pos-1), *, IOSTAT=ios) val
                IF (ios /= 0) ierr = ERR_SU2_READ
                RETURN
            END IF

            ! Skip the comma
            pos = pos + 1
        END DO

        ! Column not found in this row
        ierr = ERR_SU2_READ
    END SUBROUTINE parse_csv_value


    ! =========================================================================
    ! read_SU2_solution  --  read an SU2 ASCII restart file, convert
    !                        conservative -> primitive, extract grid coords.
    !
    ! All output arrays are ALLOCATABLE, INTENT(out).
    ! =========================================================================
    SUBROUTINE read_SU2_solution(filepath, rho, p, T, U, V, W, &
                                 Xg, Yg, Zg, n, ierr)
        CHARACTER(len=*), INTENT(IN) :: filepath
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: rho, p, T
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: U, V, W
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: Xg, Yg, Zg
        INTEGER(ik), INTENT(OUT) :: n
        INTEGER(ik), INTENT(OUT) :: ierr

        ! Locals
        INTEGER(ik) :: unit_num, ios
        CHARACTER(len=4096) :: header_line
        CHARACTER(len=4096) :: line_buf
        CHARACTER(len=64) :: col_names(MAX_SU2_COLS)
        INTEGER(ik) :: ncols
        INTEGER(ik) :: idx_x, idx_y, idx_z
        INTEGER(ik) :: idx_rho, idx_mx, idx_my, idx_mz, idx_e
        INTEGER(ik) :: nrows, i, row
        REAL(rk) :: dens, mom_x, mom_y, mom_z, energy
        REAL(rk) :: vel_u, vel_v, vel_w, ke

        ierr = 0
        n    = 0

        ! Deallocate if previously allocated
        IF (ALLOCATED(rho)) DEALLOCATE(rho)
        IF (ALLOCATED(p))   DEALLOCATE(p)
        IF (ALLOCATED(T))   DEALLOCATE(T)
        IF (ALLOCATED(U))   DEALLOCATE(U)
        IF (ALLOCATED(V))   DEALLOCATE(V)
        IF (ALLOCATED(W))   DEALLOCATE(W)
        IF (ALLOCATED(Xg))  DEALLOCATE(Xg)
        IF (ALLOCATED(Yg))  DEALLOCATE(Yg)
        IF (ALLOCATED(Zg))  DEALLOCATE(Zg)

        ! --- Open file ---
        CALL get_unit(unit_num)
        OPEN(unit=unit_num, file=TRIM(filepath), status='old', &
             action='read', iostat=ios)
        IF (ios /= 0) THEN
            ierr = ERR_SU2_OPEN
            CALL log_error(ERR_SU2_OPEN, 'File: '//TRIM(filepath))
            RETURN
        END IF

        ! --- Read and parse header line ---
        READ(unit_num, '(A)', IOSTAT=ios) header_line
        IF (ios /= 0) THEN
            ierr = ERR_SU2_HEADER
            CALL log_error(ERR_SU2_HEADER, 'File: '//TRIM(filepath))
            CLOSE(unit_num)
            RETURN
        END IF

        CALL parse_header(header_line, col_names, ncols)
        IF (ncols == 0) THEN
            ierr = ERR_SU2_HEADER
            CALL log_error(ERR_SU2_HEADER, &
                'No columns found in header: '//TRIM(filepath))
            CLOSE(unit_num)
            RETURN
        END IF

        ! --- Locate required columns by name ---
        idx_x   = find_column(col_names, ncols, 'x')
        idx_y   = find_column(col_names, ncols, 'y')
        idx_z   = find_column(col_names, ncols, 'z')
        idx_rho = find_column(col_names, ncols, 'Density')
        idx_mx  = find_column(col_names, ncols, 'Momentum_x')
        idx_my  = find_column(col_names, ncols, 'Momentum_y')
        idx_mz  = find_column(col_names, ncols, 'Momentum_z')
        idx_e   = find_column(col_names, ncols, 'Energy')

        IF (idx_x == 0 .OR. idx_y == 0 .OR. idx_z == 0) THEN
            ierr = ERR_SU2_MISSING_COLUMN
            CALL log_error(ERR_SU2_MISSING_COLUMN, &
                'Missing x/y/z in header: '//TRIM(filepath))
            CLOSE(unit_num)
            RETURN
        END IF
        IF (idx_rho == 0 .OR. idx_mx == 0 .OR. idx_my == 0 .OR. &
            idx_mz == 0 .OR. idx_e == 0) THEN
            ierr = ERR_SU2_MISSING_COLUMN
            CALL log_error(ERR_SU2_MISSING_COLUMN, &
                'Missing Density/Momentum/Energy in header: '//TRIM(filepath))
            CLOSE(unit_num)
            RETURN
        END IF

        ! --- First pass: count data rows ---
        nrows = 0
        DO
            READ(unit_num, '(A)', IOSTAT=ios) line_buf
            IF (ios /= 0) EXIT
            IF (is_data_row(line_buf)) THEN
                nrows = nrows + 1
            ELSE
                EXIT  ! hit metadata trailer
            END IF
        END DO

        IF (nrows == 0) THEN
            ierr = ERR_SU2_READ
            CALL log_error(ERR_SU2_READ, &
                'No data rows found in: '//TRIM(filepath))
            CLOSE(unit_num)
            RETURN
        END IF

        ! --- Allocate output arrays ---
        ALLOCATE(rho(nrows), p(nrows), T(nrows), &
                 U(nrows), V(nrows), W(nrows), &
                 Xg(nrows), Yg(nrows), Zg(nrows), STAT=ios)
        IF (ios /= 0) THEN
            ierr = ERR_SU2_ALLOC
            CALL log_error(ERR_SU2_ALLOC, &
                'Allocation failed for '//TRIM(INT_TO_STR(nrows))//' rows')
            CLOSE(unit_num)
            RETURN
        END IF

        ! --- Rewind and skip header, then second pass: read data ---
        REWIND(unit_num)
        READ(unit_num, '(A)', IOSTAT=ios) header_line   ! skip header

        row = 0
        DO i = 1, nrows
            READ(unit_num, '(A)', IOSTAT=ios) line_buf
            IF (ios /= 0) THEN
                ierr = ERR_SU2_READ
                CALL log_error(ERR_SU2_READ, &
                    'Unexpected EOF at row '//TRIM(INT_TO_STR(i)))
                CLOSE(unit_num)
                RETURN
            END IF

            row = row + 1

            ! Extract grid coordinates
            CALL parse_csv_value(line_buf, idx_x, Xg(row), ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing x at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF
            CALL parse_csv_value(line_buf, idx_y, Yg(row), ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing y at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF
            CALL parse_csv_value(line_buf, idx_z, Zg(row), ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing z at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF

            ! Extract conservative variables
            CALL parse_csv_value(line_buf, idx_rho, dens, ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing Density at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF
            CALL parse_csv_value(line_buf, idx_mx, mom_x, ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing Momentum_x at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF
            CALL parse_csv_value(line_buf, idx_my, mom_y, ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing Momentum_y at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF
            CALL parse_csv_value(line_buf, idx_mz, mom_z, ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing Momentum_z at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF
            CALL parse_csv_value(line_buf, idx_e, energy, ierr)
            IF (ierr /= 0) THEN
                CALL log_error(ERR_SU2_READ, &
                    'Failed parsing Energy at row '//TRIM(INT_TO_STR(row)))
                CLOSE(unit_num)
                RETURN
            END IF

            ! --- Conservative -> Primitive conversion ---
            vel_u = mom_x / dens
            vel_v = mom_y / dens
            vel_w = mom_z / dens
            ke    = 0.5_rk * dens * (vel_u*vel_u + vel_v*vel_v + vel_w*vel_w)

            rho(row) = dens
            U(row)   = vel_u
            V(row)   = vel_v
            W(row)   = vel_w
            p(row)   = (gamma_gas - 1.0_rk) * (energy - ke)
            T(row)   = p(row) / (dens * R_gas)
        END DO

        ! --- Sanity check: verify the next line is NOT another data row ---
        ! Guards against silent mesh truncation from blank lines or other
        ! anomalies that caused is_data_row to return .FALSE. mid-data.
        READ(unit_num, '(A)', IOSTAT=ios) line_buf
        IF (ios == 0) THEN
            IF (is_data_row(line_buf)) THEN
                ierr = ERR_SU2_READ
                CALL log_error(ERR_SU2_READ, &
                    'Data row found after expected end of data block '// &
                    '(row count from first pass was '// &
                    TRIM(INT_TO_STR(nrows))//'). Possible blank line '// &
                    'in data section of: '//TRIM(filepath))
                CLOSE(unit_num)
                RETURN
            END IF
        END IF

        n = nrows
        CLOSE(unit_num)

    END SUBROUTINE read_SU2_solution


    ! =========================================================================
    ! write_SU2_restart  --  convert primitive -> conservative and write an
    !                        SU2 ASCII restart file, preserving the existing
    !                        file's header and non-flow columns (PointID,
    !                        coordinates) by reading old rows in lockstep.
    !
    ! Strategy: open existing filepath for reading, open filepath.tmp for
    ! writing. Copy header verbatim, then for each data row read the old row
    ! to carry over non-flow columns and write the new conservative values in
    ! the correct columns. Copy any trailing metadata unchanged. Then rename
    ! .tmp over the original.
    ! =========================================================================
    SUBROUTINE write_SU2_restart(filepath, rho, p, T, U, V, W, n, ierr)
        CHARACTER(len=*), INTENT(IN)  :: filepath
        REAL(rk), DIMENSION(:), INTENT(IN) :: rho, p, T
        REAL(rk), DIMENSION(:), INTENT(IN) :: U, V, W
        INTEGER(ik), INTENT(IN)  :: n
        INTEGER(ik), INTENT(OUT) :: ierr

        ! Locals
        INTEGER(ik) :: unit_in, unit_out, ios
        CHARACTER(len=4096) :: header_line, line_buf, out_line
        CHARACTER(len=256)  :: temp_filepath
        CHARACTER(len=64)   :: col_names(MAX_SU2_COLS)
        CHARACTER(len=64)   :: tokens(MAX_SU2_COLS)
        INTEGER(ik) :: ncols, ntokens
        INTEGER(ik) :: idx_rho, idx_mx, idx_my, idx_mz, idx_e
        INTEGER(ik) :: i, j
        REAL(rk) :: dens, mom_x, mom_y, mom_z, energy
        CHARACTER(len=32) :: val_str

        ierr = 0
        temp_filepath = TRIM(filepath)//'.tmp'

        ! --- Open existing file for reading ---
        CALL get_unit(unit_in)
        OPEN(unit=unit_in, file=TRIM(filepath), status='old', &
             action='read', iostat=ios)
        IF (ios /= 0) THEN
            ierr = ERR_SU2_OPEN
            CALL log_error(ERR_SU2_OPEN, 'File: '//TRIM(filepath))
            RETURN
        END IF

        ! --- Open temp file for writing ---
        CALL get_unit(unit_out)
        OPEN(unit=unit_out, file=TRIM(temp_filepath), status='replace', &
             action='write', iostat=ios)
        IF (ios /= 0) THEN
            ierr = ERR_SU2_OPEN
            CALL log_error(ERR_SU2_OPEN, 'File: '//TRIM(temp_filepath))
            CLOSE(unit_in)
            RETURN
        END IF

        ! --- Read and write header verbatim ---
        READ(unit_in, '(A)', IOSTAT=ios) header_line
        IF (ios /= 0) THEN
            ierr = ERR_SU2_HEADER
            CALL log_error(ERR_SU2_HEADER, 'File: '//TRIM(filepath))
            CLOSE(unit_in); CLOSE(unit_out)
            RETURN
        END IF
        WRITE(unit_out, '(A)', IOSTAT=ios) TRIM(header_line)
        IF (ios /= 0) THEN
            ierr = ERR_SU2_WRITE
            CALL log_error(ERR_SU2_WRITE, 'Header write failed')
            CLOSE(unit_in); CLOSE(unit_out)
            RETURN
        END IF

        ! --- Parse header to find conservative column indices ---
        CALL parse_header(header_line, col_names, ncols)
        idx_rho = find_column(col_names, ncols, 'Density')
        idx_mx  = find_column(col_names, ncols, 'Momentum_x')
        idx_my  = find_column(col_names, ncols, 'Momentum_y')
        idx_mz  = find_column(col_names, ncols, 'Momentum_z')
        idx_e   = find_column(col_names, ncols, 'Energy')

        IF (idx_rho == 0 .OR. idx_mx == 0 .OR. idx_my == 0 .OR. &
            idx_mz == 0 .OR. idx_e == 0) THEN
            ierr = ERR_SU2_MISSING_COLUMN
            CALL log_error(ERR_SU2_MISSING_COLUMN, &
                'Missing flow columns in header for write: '//TRIM(filepath))
            CLOSE(unit_in); CLOSE(unit_out)
            RETURN
        END IF

        ! --- Process data rows: read old, replace flow columns, write new ---
        DO i = 1, n
            READ(unit_in, '(A)', IOSTAT=ios) line_buf
            IF (ios /= 0) THEN
                ierr = ERR_SU2_READ
                CALL log_error(ERR_SU2_READ, &
                    'Unexpected EOF at row '//TRIM(INT_TO_STR(i)))
                CLOSE(unit_in); CLOSE(unit_out)
                RETURN
            END IF

            ! Tokenize the old row
            CALL tokenize_csv(line_buf, tokens, ntokens)

            ! Compute conservative values from primitive
            dens  = rho(i)
            mom_x = rho(i) * U(i)
            mom_y = rho(i) * V(i)
            mom_z = rho(i) * W(i)
            energy = p(i) / (gamma_gas - 1.0_rk) + &
                     0.5_rk * rho(i) * (U(i)*U(i) + V(i)*V(i) + W(i)*W(i))

            ! Replace flow columns in the token array
            WRITE(val_str, '(ES23.15E3)') dens
            tokens(idx_rho) = ADJUSTL(val_str)
            WRITE(val_str, '(ES23.15E3)') mom_x
            tokens(idx_mx)  = ADJUSTL(val_str)
            WRITE(val_str, '(ES23.15E3)') mom_y
            tokens(idx_my)  = ADJUSTL(val_str)
            WRITE(val_str, '(ES23.15E3)') mom_z
            tokens(idx_mz)  = ADJUSTL(val_str)
            WRITE(val_str, '(ES23.15E3)') energy
            tokens(idx_e)   = ADJUSTL(val_str)

            ! Reassemble the line
            out_line = ''
            DO j = 1, ntokens
                IF (j == 1) THEN
                    out_line = TRIM(tokens(j))
                ELSE
                    out_line = TRIM(out_line)//','//TRIM(tokens(j))
                END IF
            END DO

            WRITE(unit_out, '(A)', IOSTAT=ios) TRIM(out_line)
            IF (ios /= 0) THEN
                ierr = ERR_SU2_WRITE
                CALL log_error(ERR_SU2_WRITE, &
                    'Write failed at row '//TRIM(INT_TO_STR(i)))
                CLOSE(unit_in); CLOSE(unit_out)
                RETURN
            END IF
        END DO

        ! --- Copy any trailing metadata lines unchanged ---
        DO
            READ(unit_in, '(A)', IOSTAT=ios) line_buf
            IF (ios /= 0) EXIT
            WRITE(unit_out, '(A)', IOSTAT=ios) TRIM(line_buf)
            IF (ios /= 0) THEN
                ierr = ERR_SU2_WRITE
                CALL log_error(ERR_SU2_WRITE, 'Metadata copy failed')
                CLOSE(unit_in); CLOSE(unit_out)
                RETURN
            END IF
        END DO

        CLOSE(unit_in)
        CLOSE(unit_out)

        ! --- Atomic rename: .tmp over original ---
        BLOCK
            INTEGER(C_INT) :: rc
            rc = c_rename(TRIM(temp_filepath)//C_NULL_CHAR, &
                          TRIM(filepath)//C_NULL_CHAR)
            IF (rc /= 0_C_INT) THEN
                ierr = ERR_SU2_RENAME
                CALL log_error(ERR_SU2_RENAME, &
                    'rename() failed: '//TRIM(temp_filepath)// &
                    ' -> '//TRIM(filepath))
                RETURN
            END IF
        END BLOCK

    END SUBROUTINE write_SU2_restart


    ! -------------------------------------------------------------------------
    ! tokenize_csv  --  split a comma-separated line into token strings.
    ! -------------------------------------------------------------------------
    SUBROUTINE tokenize_csv(line, tokens, ntokens)
        CHARACTER(len=*), INTENT(IN) :: line
        CHARACTER(len=64), INTENT(OUT) :: tokens(MAX_SU2_COLS)
        INTEGER(ik), INTENT(OUT) :: ntokens

        INTEGER :: pos, start, llen

        ntokens = 0
        llen    = LEN_TRIM(line)
        pos     = 1

        DO WHILE (pos <= llen)
            start = pos
            ! Advance to next comma or end
            DO WHILE (pos <= llen .AND. line(pos:pos) /= ',')
                pos = pos + 1
            END DO

            ntokens = ntokens + 1
            IF (ntokens <= MAX_SU2_COLS) THEN
                IF (start <= pos - 1) THEN
                    tokens(ntokens) = ADJUSTL(line(start:pos-1))
                ELSE
                    tokens(ntokens) = ''   ! empty field (trailing comma)
                END IF
            END IF

            ! Skip the comma
            pos = pos + 1
        END DO

        ! Handle trailing comma: if the line ends with ',' the loop above
        ! exits with pos = llen+2, but the empty field after the final
        ! comma was never captured. Detect and add it.
        IF (llen > 0 .AND. line(llen:llen) == ',') THEN
            ntokens = ntokens + 1
            IF (ntokens <= MAX_SU2_COLS) tokens(ntokens) = ''
        END IF
    END SUBROUTINE tokenize_csv

END MODULE SU2_IO
