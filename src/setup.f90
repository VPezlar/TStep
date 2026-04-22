! =============================================================================
! setup  --  reads inputs/inputs.in and provides small utility helpers.
!
! `configurationRead` parses three Fortran namelists into the globals in
! `variables`:
!   &General  -- solver-agnostic settings (flow_format, output_file, eps_0, ...)
!   &Arnoldi  -- eigenvalue-solver parameters
!   &<solver> -- format-specific block, chosen via flow_format (e.g. OpenFOAM)
!
! Path pipeline (flat, no module-globals for derived paths):
!   validate_paths(ierr)        -- sanity-check the three raw inputs from the
!                                  &OpenFOAM namelist (non-empty grid path,
!                                  trailing slash on both directories).
!   seed_stability_initial(ierr)-- mkdir -p <stability_dir>/1/ and cp -a the
!                                  four baseflow fields (p, rho, T, U) from
!                                  <baseflow_field> into it. Required before
!                                  the first write_flowfield call because
!                                  write_OF_scalars/vectors preserve the
!                                  OpenFOAM ASCII header of a PRE-EXISTING
!                                  target file and only rewrite the body.
!   promote_to_initial_state(ierr) -- cp -af the four OpenFOAM fields from
!                                  <stability_dir>/<1+TTime>/ back over
!                                  <stability_dir>/1/. Used once after the
!                                  CFD warmup pass to turn `1/` into
!                                  `q0 = F(q_raw)`, the discretely-consistent
!                                  base state used by the Arnoldi loop.
!   stability_time_dir(t)       -- '<stability_dir>/<time_to_str(t)>/';
!                                  callers build <stab>/1/ or <stab>/<1+TTime>/
!                                  on demand without any stored globals.
!   stability_output_path(name) -- '<stability_dir>/<name>'; used by CSV/
!                                  eigendata writers to land all TStep
!                                  artefacts under the stability case root
!                                  instead of scattering them around disk.
!
! Other helpers:
!   get_unit(u)      -- return an unused Fortran unit number in [10, 99]
!   INT_TO_STR(val)  -- cheap integer-to-string for error messages
!   time_to_str(t)   -- format a real as OpenFOAM 'general, precision 6' would
! =============================================================================
MODULE setup
    USE accuracy
    USE variables
    USE error_handling

    IMPLICIT NONE

    PUBLIC :: configurationRead, validate_paths, seed_stability_initial, &
              promote_to_initial_state, stability_time_dir, &
              stability_output_path, clear_endpoint_folder, &
              find_endpoint_folder, get_unit, INT_TO_STR, time_to_str

CONTAINS

    SUBROUTINE configurationRead(ierr)
        INTEGER(ik), INTENT(out) :: ierr
        INTEGER(ik) :: unit_num, status_id
        CHARACTER(len=400) :: message

        ! General namelist (solver-agnostic)
        NAMELIST / General / flow_format, &
                             eps_0, &
                             COMMAND_RUN

        ! Arnoldi-specific namelist
        NAMELIST / Arnoldi / krylov_size, &
                             frechet_order, &
                             eigenvalue_sort_by, &
                             TTime, &
                             num_threads, &
                             N_eig_write

        ! OpenFOAM-specific namelist (raw path inputs)
        NAMELIST / OpenFOAM / N_HEADER_grid, &
                              N_HEADER_var, &
                              baseflow_grid, &
                              baseflow_field, &
                              stability_dir

        ierr = 0
        CALL get_unit(unit_num)
        OPEN(unit_num, file="../inputs/inputs.in", status="old", &
             iostat=status_id, iomsg=message)

        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_FILE_OPEN
            CALL log_error(ERR_SETUP_FILE_OPEN, &
                'File: ../inputs/inputs.in - '//TRIM(message))
            RETURN
        END IF

        ! Read General namelist
        READ(unit_num, nml=General, iostat=status_id, iomsg=message)
        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_NAMELIST_READ
            CALL log_error(ERR_SETUP_NAMELIST_READ, &
                'General namelist - '//TRIM(message))
            CLOSE(unit_num)
            RETURN
        END IF

        ! Read Arnoldi namelist
        READ(unit_num, nml=Arnoldi, iostat=status_id, iomsg=message)
        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_NAMELIST_READ
            CALL log_error(ERR_SETUP_NAMELIST_READ, &
                'Arnoldi namelist - '//TRIM(message))
            CLOSE(unit_num)
            RETURN
        END IF

        ! Read solver-specific namelist based on flow_format
        IF (TRIM(flow_format) == 'OpenFOAM') THEN
            READ(unit_num, nml=OpenFOAM, iostat=status_id, iomsg=message)
            IF (status_id /= 0) THEN
                ierr = ERR_SETUP_NAMELIST_READ
                CALL log_error(ERR_SETUP_NAMELIST_READ, &
                    'OpenFOAM namelist - '//TRIM(message))
                CLOSE(unit_num)
                RETURN
            END IF
        ELSE
            ! Future solvers can be added here
            ierr = ERR_SETUP_INVALID_PARAM
            CALL log_error(ERR_SETUP_INVALID_PARAM, &
                'Unsupported flow_format: '//TRIM(flow_format))
            CLOSE(unit_num)
            RETURN
        END IF

        CLOSE(unit_num)

    END SUBROUTINE configurationRead


    ! -------------------------------------------------------------------------
    ! validate_paths  --  sanity-check the three raw &OpenFOAM path inputs.
    !
    ! Nothing is derived or stashed in module globals. Call sites compose
    ! time-folder paths on demand via stability_time_dir(t). This keeps the
    ! path story obvious: what you set in inputs.in is exactly what TStep
    ! reads/writes, give or take the '1/' or '<1+TTime>/' suffix.
    !
    ! Fails with ERR_SETUP_INVALID_PARAM when:
    !   - baseflow_grid is empty (it is a FILE path, the grid cell-centers).
    !   - baseflow_field or stability_dir is missing a trailing '/' (required
    !     so concatenation yields a valid directory path).
    ! -------------------------------------------------------------------------
    SUBROUTINE validate_paths(ierr)
        INTEGER(ik), INTENT(OUT) :: ierr
        INTEGER :: ng, nf, ns

        ierr = 0

        ng = LEN_TRIM(baseflow_grid)
        IF (ng == 0) THEN
            ierr = ERR_SETUP_INVALID_PARAM
            CALL log_error(ERR_SETUP_INVALID_PARAM, &
                'baseflow_grid must be a non-empty path to the grid C file')
            RETURN
        END IF

        nf = LEN_TRIM(baseflow_field)
        IF (nf == 0 .OR. baseflow_field(nf:nf) /= '/') THEN
            ierr = ERR_SETUP_INVALID_PARAM
            CALL log_error(ERR_SETUP_INVALID_PARAM, &
                'baseflow_field must end with "/": '//TRIM(baseflow_field))
            RETURN
        END IF

        ns = LEN_TRIM(stability_dir)
        IF (ns == 0 .OR. stability_dir(ns:ns) /= '/') THEN
            ierr = ERR_SETUP_INVALID_PARAM
            CALL log_error(ERR_SETUP_INVALID_PARAM, &
                'stability_dir must end with "/": '//TRIM(stability_dir))
            RETURN
        END IF

        ! Ensure the standard output subfolder exists under stability_dir.
        ! TStep deposits flowfield.csv / eigenvalues.dat / eigenvectors.dat
        ! in <stability_dir>/output/, composed internally (no user knob).
        CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(stability_dir)//'output/', &
                                  wait=.TRUE.)

        WRITE(*,'(A)')         '---------------------------------------'
        WRITE(*,'(A)')         'TStep validated paths:'
        WRITE(*,'(A,A)')       '  baseflow_grid  = ', TRIM(baseflow_grid)
        WRITE(*,'(A,A)')       '  baseflow_field = ', TRIM(baseflow_field)
        WRITE(*,'(A,A)')       '  stability_dir  = ', TRIM(stability_dir)
        WRITE(*,'(A,A)')       '  initial seed   -> ', &
            TRIM(stability_time_dir(TSTEP_INITIAL_TIME))
        WRITE(*,'(A,A)')       '  solver result  <- ', &
            TRIM(stability_time_dir(TSTEP_INITIAL_TIME + TTime))
        WRITE(*,'(A,F12.6,A)') '  TSTEP_INITIAL_TIME= ', TSTEP_INITIAL_TIME, &
                               '   (controlDict.startTime must match)'
        WRITE(*,'(A,F12.6,A)') '  end time (1+TTime)= ', &
                               TSTEP_INITIAL_TIME+TTime, &
                               '   (controlDict.endTime must match)'
        WRITE(*,'(A)')         '---------------------------------------'
    END SUBROUTINE validate_paths


    ! -------------------------------------------------------------------------
    ! seed_stability_initial  --  prepare <stability_dir>/1/ for the first
    !                              write_flowfield call.
    !
    ! Why this exists: write_OF_scalars and write_OF_vectors in OpenFOAM_IO
    ! open the target with status='old' -- they preserve the OpenFOAM ASCII
    ! header and only rewrite the numeric body. So <stability_dir>/1/{p,rho,
    ! T,U} must exist on disk before TStep attempts the write.
    !
    ! Procedure (shell-assisted via EXECUTE_COMMAND_LINE):
    !   1. mkdir -p <stability_dir>/1/
    !   2. cp -a <baseflow_field>{p,rho,T,U} <stability_dir>/1/
    !
    ! Net effect: the '1/' folder now holds full OpenFOAM files whose header
    ! matches the baseflow (correct dimensions, boundary fields, etc.) and
    ! whose body is about to be overwritten with (base + perturbation).
    ! -------------------------------------------------------------------------
    SUBROUTINE seed_stability_initial(ierr)
        INTEGER(ik), INTENT(OUT) :: ierr
        CHARACTER(len=512) :: cmd
        CHARACTER(len=256) :: init_dir
        INTEGER :: cmdstat, exitstat

        ierr = 0
        init_dir = TRIM(stability_time_dir(TSTEP_INITIAL_TIME))

        ! 1. Create the destination time folder (idempotent).
        cmd = 'mkdir -p ' // TRIM(init_dir)
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE., &
                                  exitstat=exitstat, cmdstat=cmdstat)
        IF (cmdstat /= 0 .OR. exitstat /= 0) THEN
            ierr = ERR_SETUP_FILE_OPEN
            CALL log_error(ERR_SETUP_FILE_OPEN, &
                'mkdir -p failed: '//TRIM(init_dir))
            RETURN
        END IF

        ! 2. Seed p, rho, T, U with the baseflow templates. `cp -a` preserves
        !    permissions/timestamps; -f overwrites any stale files from a
        !    previous TStep run so re-runs are deterministic.
        cmd = 'cp -af ' // TRIM(baseflow_field) // 'p '    // &
                          TRIM(baseflow_field) // 'rho '  // &
                          TRIM(baseflow_field) // 'T '    // &
                          TRIM(baseflow_field) // 'U '    // &
                          TRIM(init_dir)
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE., &
                                  exitstat=exitstat, cmdstat=cmdstat)
        IF (cmdstat /= 0 .OR. exitstat /= 0) THEN
            ierr = ERR_SETUP_FILE_OPEN
            CALL log_error(ERR_SETUP_FILE_OPEN, &
                'cp -af failed: '//TRIM(baseflow_field)//'{p,rho,T,U} -> '//TRIM(init_dir))
            RETURN
        END IF

        WRITE(*,'(A,A)') 'Seeded initial state: ', TRIM(init_dir)
    END SUBROUTINE seed_stability_initial


    ! -------------------------------------------------------------------------
    ! stability_time_dir  --  compose '<stability_dir>/<time_to_str(t)>/'.
    !
    ! Callers use this for the two time folders TStep touches:
    !   stability_time_dir(TSTEP_INITIAL_TIME)        -- '<stab>/1/'
    !   stability_time_dir(TSTEP_INITIAL_TIME+TTime)  -- '<stab>/1.<..>/'
    ! -------------------------------------------------------------------------
    PURE FUNCTION stability_time_dir(t) RESULT(d)
        REAL(rk), INTENT(IN) :: t
        CHARACTER(len=256)   :: d
        d = TRIM(stability_dir) // TRIM(time_to_str(t)) // '/'
    END FUNCTION stability_time_dir


    ! -------------------------------------------------------------------------
    ! promote_to_initial_state  --  make the CFD solver's endpoint output the
    !                                new '1/' snapshot.
    !
    ! Copies the four OpenFOAM field files (p, rho, T, U) from
    ! <stability_dir>/<1+TTime>/ back over <stability_dir>/1/. Same shell
    ! trick as seed_stability_initial: a single `cp -af` list so the headers
    ! and numeric bodies arrive together.
    !
    ! Used once during init, right after the CFD warmup pass, so that `1/`
    ! holds q0 = F(q_raw) (the discretely-consistent base state) for the
    ! subsequent Arnoldi loop.
    ! -------------------------------------------------------------------------
    SUBROUTINE promote_to_initial_state(ierr)
        INTEGER(ik), INTENT(OUT) :: ierr
        CHARACTER(len=512) :: cmd
        CHARACTER(len=256) :: init_dir, end_dir
        INTEGER :: cmdstat, exitstat

        ierr = 0
        init_dir = TRIM(stability_time_dir(TSTEP_INITIAL_TIME))
        CALL find_endpoint_folder(end_dir, ierr)
        IF (ierr /= 0) RETURN

        cmd = 'cp -af ' // TRIM(end_dir) // 'p '    // &
                          TRIM(end_dir) // 'rho '  // &
                          TRIM(end_dir) // 'T '    // &
                          TRIM(end_dir) // 'U '    // &
                          TRIM(init_dir)
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE., &
                                  exitstat=exitstat, cmdstat=cmdstat)
        IF (cmdstat /= 0 .OR. exitstat /= 0) THEN
            ierr = ERR_SETUP_FILE_OPEN
            CALL log_error(ERR_SETUP_FILE_OPEN, &
                'cp -af failed: '//TRIM(end_dir)//'{p,rho,T,U} -> '//TRIM(init_dir))
            RETURN
        END IF

        WRITE(*,'(A,A,A,A)') 'Promoted solver output ', TRIM(end_dir), &
                             ' -> ', TRIM(init_dir)
    END SUBROUTINE promote_to_initial_state


    ! -------------------------------------------------------------------------
    ! stability_output_path  --  compose '<stability_dir>/output/<filename>'.
    ! The 'output/' subfolder is a standardised location composed
    ! internally; users configure only stability_dir. validate_paths
    ! mkdir -p's it on startup so writes never fail on a missing parent.
    !
    ! Used by write_output / write_eigendata so every TStep artefact lands
    ! under the stability case root (sacred baseflow folders are never
    ! touched). `stability_dir` is guaranteed to exist by validate_paths.
    ! -------------------------------------------------------------------------
    PURE FUNCTION stability_output_path(filename) RESULT(p)
        CHARACTER(len=*), INTENT(IN) :: filename
        CHARACTER(len=256)           :: p
        p = TRIM(stability_dir) // 'output/' // TRIM(filename)
    END FUNCTION stability_output_path

    ! -------------------------------------------------------------------------
    ! clear_endpoint_folder  --  wipe every non-initial numeric time folder
    ! under stability_dir so the solver is guaranteed a clean slate for its
    ! endpoint write.
    !
    ! Rationale: OpenFOAM's time loop accumulates IEEE-754 round-off. For
    ! TTime = 0.1 with deltaT = 1e-4, after ~1000 steps the internal time is
    ! 1.099999999999989 (not 1.1 exactly). Depending on timePrecision and
    ! pre-existing folders, OpenFOAM may write '<stab>/1.1/' OR
    ! '<stab>/1.099999999999989/' -- and can even auto-bump timePrecision to
    ! disambiguate. TStep is not in the business of predicting that name
    ! anymore: we clear all candidate endpoint folders first, then pick
    ! whatever single folder appears after the solver call via
    ! find_endpoint_folder() below.
    !
    ! The glob deliberately excludes '<stab>/1/' (the initial state) and
    ! non-numeric folders (constant, system, output, ...). Safe because the
    ! deletion targets are strictly inside stability_dir.
    ! -------------------------------------------------------------------------
    SUBROUTINE clear_endpoint_folder(ierr)
        INTEGER(ik), INTENT(OUT) :: ierr
        CHARACTER(len=1024) :: cmd
        ierr = 0
        ! find <stab> -maxdepth 1 -type d -regex '<stab>/[0-9][0-9.]*' -not
        ! -path '<stab>/1' -exec rm -rf {} +
        cmd = 'find '//TRIM(stability_dir)//&
              ' -maxdepth 1 -type d -regextype posix-extended'//&
              ' -regex '''//TRIM(stability_dir)//'[0-9]+(\.[0-9]+)?'''//&
              ' -not -path '''//TRIM(stability_dir)//'1'''//&
              ' -exec rm -rf {} +'
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE.)
    END SUBROUTINE clear_endpoint_folder


    ! -------------------------------------------------------------------------
    ! find_endpoint_folder  --  return the latest numeric-named subfolder of
    ! stability_dir, excluding the initial '1/' seed.
    !
    ! Used as the source path for every endpoint read. Fully decouples TStep
    ! from OpenFOAM's timePrecision / IEEE-round-off dance: whatever the
    ! solver decided to call its output folder ('1.1/', '1.099999999999989/',
    ! etc.), we pick it up. Returns an empty string + ierr on failure so the
    ! caller can surface a clear error.
    !
    ! Implementation: shell pipeline writes the folder name to a tmp file,
    ! Fortran reads it back. ls+awk+sort -g+tail -n1 gives us the largest
    ! numeric folder name > '1/'.
    ! -------------------------------------------------------------------------
    SUBROUTINE find_endpoint_folder(path, ierr)
        CHARACTER(len=256), INTENT(OUT) :: path
        INTEGER(ik),        INTENT(OUT) :: ierr
        CHARACTER(len=1024) :: cmd
        CHARACTER(len=256)  :: folder_name
        CHARACTER(len=256)  :: tmpfile
        INTEGER             :: u, io_stat

        ierr        = 0
        path        = ''
        folder_name = ''
        tmpfile     = '/tmp/tstep_endpoint.txt'

        ! List <stab>, keep numeric names != '1', sort by value, pick max.
        ! awk pattern: /^[0-9]/ matches numeric-starting names; $0 != "1"
        ! excludes the initial-state folder.
        cmd = 'cd '//TRIM(stability_dir)//' && ls -1 2>/dev/null'//&
              ' | awk ''/^[0-9]/ && $0 != "1"'''//&
              ' | sort -g | tail -n 1 > '//TRIM(tmpfile)
        CALL EXECUTE_COMMAND_LINE(TRIM(cmd), wait=.TRUE.)

        CALL get_unit(u)
        OPEN(u, file=TRIM(tmpfile), status='old', action='read', iostat=io_stat)
        IF (io_stat == 0) THEN
            READ(u, '(A)', iostat=io_stat) folder_name
            CLOSE(u)
        END IF

        IF (LEN_TRIM(folder_name) == 0) THEN
            ierr = ERR_SETUP_INVALID_PARAM
            CALL log_error(ERR_SETUP_INVALID_PARAM, &
                'find_endpoint_folder: no post-solver time folder found under '&
                //TRIM(stability_dir)//' (solver did not write, or controlDict wrong)')
            RETURN
        END IF

        path = TRIM(stability_dir)//TRIM(folder_name)//'/'
    END SUBROUTINE find_endpoint_folder


    ! -------------------------------------------------------------------------
    ! time_to_str  --  format a real number the way OpenFOAM does under
    !                  'timeFormat general, timePrecision 15'.
    !
    ! Bulletproof choice: 15 fractional digits sits right at the IEEE 754
    ! double-precision limit (REAL64 has ~15.95 sig figs), so Fortran's
    ! round-half-to-even absorbs the trailing-bit noise that IEEE addition
    ! introduces (e.g. 1.0d0 + 0.1d0 = 1.10000000000000008881...). After
    ! the trailing-zero strip below, such values collapse back to their
    ! canonical form ('1.1'), which is exactly what OpenFOAM writes at
    ! timePrecision 15 in 'general' format.
    !
    ! Rule:
    !   - For |t| >= 1e-4 or t == 0:  fixed-point, up to 15 decimal digits,
    !     trailing zeros (and the decimal point if bare) stripped.
    !     Examples: 1.0 -> '1',  1.1 -> '1.1',  0.02 -> '0.02',
    !               1.000001 -> '1.000001',
    !               1.123456789012345 -> '1.123456789012345'.
    !   - For 0 < |t| < 1e-4:  scientific. NOT supported by TStep today;
    !     we abort via empty-string return. Callers must stay in range.
    !
    ! Matching OpenFOAM's convention exactly lets us predict the solver-output
    ! folder name without running OpenFOAM first. The controlDict MUST use
    ! timePrecision 15 (recommended) or any lower value for which the
    ! endpoint time has no digits that would be truncated.
    ! -------------------------------------------------------------------------
    PURE FUNCTION time_to_str(t) RESULT(s)
        REAL(rk), INTENT(IN) :: t
        CHARACTER(len=32)    :: s
        CHARACTER(len=32)    :: buf
        INTEGER              :: i, dot

        IF (ABS(t) >= 1.0e-4_rk .OR. t == 0.0_rk) THEN

            ! F0.6 can produce '.02' (no leading zero) on some compilers.
            ! OpenFOAM always writes '0.02', so reinsert the leading zero.
            WRITE(buf, '(F0.15)') t
            IF (buf(1:1) == '.') THEN
                buf = '0' // buf(1:LEN_TRIM(buf))
            ELSE IF (buf(1:2) == '-.') THEN
                buf = '-0' // buf(2:LEN_TRIM(buf))
            END IF

            ! Strip trailing zeros after the decimal point; then strip a
            ! bare trailing '.' if it was left exposed (e.g. '1.000000' -> '1').
            dot = INDEX(buf, '.')
            IF (dot > 0) THEN
                i = LEN_TRIM(buf)
                DO WHILE (i > dot .AND. buf(i:i) == '0')
                    i = i - 1
                END DO
                IF (buf(i:i) == '.') i = i - 1
                s = buf(1:i)
            ELSE
                s = TRIM(buf)
            END IF

        ELSE
            ! Out of supported range. Return empty string; calling code that
            ! builds paths from the result will fail loudly at path use.
            s = ''
        END IF
    END FUNCTION time_to_str


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


    ! Helper function to convert integer to string
    FUNCTION INT_TO_STR(val) RESULT(str)
        INTEGER(ik), INTENT(IN) :: val
        CHARACTER(len=20) :: str
        WRITE(str, '(I0)') val
    END FUNCTION INT_TO_STR

END MODULE setup