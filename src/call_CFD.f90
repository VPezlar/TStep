! =============================================================================
! call_CFD  --  thin wrapper around EXECUTE_COMMAND_LINE for invoking the
!               external CFD solver (e.g. `rhoCentralFoam -case <dir>`).
!
! `run_simulation(COMMAND_STRING, ERROR_STATUS)` blocks until the child
! process finishes, then maps its exit status onto TStep's error codes:
!   - failure to launch     -> ERR_CMD_LAUNCH_FAILED  (CMDSTAT /= 0)
!   - non-zero solver exit  -> ERR_CMD_NONZERO_EXIT   (EXITSTAT /= 0)
! The command string is whatever the user set in inputs.in under COMMAND_RUN.
!
! This module is the only place TStep hands control to an external process;
! it is the future hook for the Frechet-derivative matvec in Arnoldi.f90.
! =============================================================================
MODULE call_CFD
    USE accuracy
    USE error_handling
    USE setup, ONLY: INT_TO_STR
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: run_simulation

CONTAINS

    SUBROUTINE run_simulation(COMMAND_STRING, ERROR_STATUS)
        ! Arguments
        CHARACTER(LEN=*), INTENT(IN) :: COMMAND_STRING
        INTEGER(ik), INTENT(OUT)         :: ERROR_STATUS

        ! Local Variables (using default kind for integers since they are OS-specific)
        INTEGER(ik) :: EXIT_STAT_VAL    ! Stores the exit code of the external program (0 for success)
        INTEGER(ik) :: CMD_STAT_VAL     ! Stores the status of the command launch itself (0 if successful)

        ! --- Initialization ---
        ERROR_STATUS = 0

        WRITE(*,*) '---------------------------------------'
        WRITE(*,*) 'ATTEMPTING TO RUN EXTERNAL COMMAND:'
        WRITE(*,*) 'COMMAND: ', TRIM(COMMAND_STRING)
        
        ! --- Execution ---
        ! EXECUTE_COMMAND_LINE pauses the Fortran code until the external command finishes (WAIT=.TRUE. is default)
        CALL EXECUTE_COMMAND_LINE( &
            COMMAND    = COMMAND_STRING, &
            EXITSTAT   = EXIT_STAT_VAL,  &
            CMDSTAT    = CMD_STAT_VAL    &
        )

        ! --- Check Status ---
        IF (CMD_STAT_VAL /= 0) THEN
            ! Command failed to launch (e.g., 'mpirun' not found)
            ERROR_STATUS = ERR_CMD_LAUNCH_FAILED
            CALL log_error(ERR_CMD_LAUNCH_FAILED, 'CMDSTAT: '//TRIM(ADJUSTL(INT_TO_STR(CMD_STAT_VAL))))
            RETURN
        END IF
        
        IF (EXIT_STAT_VAL /= 0) THEN
            ! Command launched but finished with a non-zero (error) code
            ERROR_STATUS = ERR_CMD_NONZERO_EXIT
            CALL log_error(ERR_CMD_NONZERO_EXIT, 'Exit code: '//TRIM(ADJUSTL(INT_TO_STR(EXIT_STAT_VAL))))
            RETURN
        END IF
        
        ! Success
        WRITE(*,*) 'SUCCESS: Command finished successfully (Exit Code 0).'
        
    END SUBROUTINE run_simulation

END MODULE call_CFD
