PROGRAM main
    USE accuracy
    USE setup
    USE read_flow
    USE write_output
    USE call_CFD

    IMPLICIT NONE

    REAL(rk), DIMENSION(:), ALLOCATABLE :: rho_in, p_in, T_in, U_in, V_in, W_in, Xgrid, Ygrid, Zgrid
    INTEGER(ik) :: data_count
    INTEGER(ik) :: error_status, STATUS_CODE

        ! Read configuration
    CALL configurationRead(error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'FATAL: Configuration read failed.'
        RETURN
    END IF

        ! Call the new subroutine
    CALL run_simulation(COMMAND_RUN, STATUS_CODE)

    IF (STATUS_CODE /= 0) THEN
        WRITE(*,*) 'EXECUTION STOPPED DUE TO EXTERNAL PROGRAM FAILURE.'
        STOP
    END IF

    ! Continue with other code if STATUS_CODE is 0
    WRITE(*,*) 'PROCEEDING TO NEXT STEP.'

    ! --- Main Execution ---
    ! Read flowfield data
    CALL read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, &
                        Xgrid, Ygrid, Zgrid, data_count, error_status)
    
    IF (error_status /= 0) THEN
        WRITE(*,*) 'FATAL: Flowfield read failed. Exiting.'
        CALL cleanup_allocations()
        STOP 1
    END IF

    ! Write flowfield data to file
    CALL write_flowfield_data(Xgrid, Ygrid, Zgrid, rho_in, p_in, T_in, &
                              U_in, V_in, W_in, data_count, error_status)
    IF (error_status /= 0) THEN
        WRITE(*,*) 'FATAL: Failed to write flowfield data. Exiting.'
        CALL cleanup_allocations()
        STOP 2
    END IF


    ! Cleanup allocated memory
    CALL cleanup_allocations()

CONTAINS

    SUBROUTINE cleanup_allocations()
        IF (ALLOCATED(p_in)) DEALLOCATE(p_in)
        IF (ALLOCATED(rho_in)) DEALLOCATE(rho_in)
        IF (ALLOCATED(T_in)) DEALLOCATE(T_in)
        IF (ALLOCATED(U_in)) DEALLOCATE(U_in)
        IF (ALLOCATED(V_in)) DEALLOCATE(V_in)
        IF (ALLOCATED(W_in)) DEALLOCATE(W_in)
        IF (ALLOCATED(Xgrid)) DEALLOCATE(Xgrid)
        IF (ALLOCATED(Ygrid)) DEALLOCATE(Ygrid)
        IF (ALLOCATED(Zgrid)) DEALLOCATE(Zgrid)
    END SUBROUTINE cleanup_allocations

END PROGRAM main
