MODULE call_CFD
    USE accuracy
    USE error_handling
    
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
    
    ! Helper function to convert integer to string
    FUNCTION INT_TO_STR(val) RESULT(str)
        INTEGER(ik), INTENT(IN) :: val
        CHARACTER(len=20) :: str
        WRITE(str, '(I0)') val
    END FUNCTION INT_TO_STR

END MODULE call_CFD
