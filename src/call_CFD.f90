MODULE call_CFD
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: run_simulation

CONTAINS

    SUBROUTINE run_simulation(COMMAND_STRING, ERROR_STATUS)
        ! Arguments
        CHARACTER(LEN=*), INTENT(IN) :: COMMAND_STRING
        INTEGER, INTENT(OUT)         :: ERROR_STATUS

        ! Local Variables (using default kind for integers since they are OS-specific)
        INTEGER :: EXIT_STAT_VAL    ! Stores the exit code of the external program (0 for success)
        INTEGER :: CMD_STAT_VAL     ! Stores the status of the command launch itself (0 if successful)

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
            WRITE(*,*) 'FATAL: OS failed to launch command. CMDSTAT:', CMD_STAT_VAL
            ERROR_STATUS = 1 ! Use a unique error code for launch failure
            RETURN
        END IF
        
        IF (EXIT_STAT_VAL /= 0) THEN
            ! Command launched but finished with a non-zero (error) code
            WRITE(*,*) 'FATAL: Command finished with non-zero exit code:', EXIT_STAT_VAL
            ERROR_STATUS = 2 ! Use a unique error code for execution failure
            RETURN
        END IF
        
        ! Success
        WRITE(*,*) 'SUCCESS: Command finished successfully (Exit Code 0).'
        
    END SUBROUTINE run_simulation

END MODULE call_CFD