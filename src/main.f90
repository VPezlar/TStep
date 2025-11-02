PROGRAM main
    USE accuracy
    USE setup
    USE read_flow
    USE write_output
    USE call_CFD
    USE random_disturbance
    USE error_handling

    IMPLICIT NONE

    REAL(rk), DIMENSION(:), ALLOCATABLE :: rho_in, p_in, T_in, U_in, V_in, W_in, Xgrid, Ygrid, Zgrid, pert_0
    INTEGER(ik) :: data_count
    INTEGER(ik) :: error_status, STATUS_CODE

    ! Read configuration
    CALL configurationRead(error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_CONFIG)
        STOP ERR_MAIN_CONFIG
    END IF

    ! Execute external command
    CALL run_simulation(COMMAND_RUN, STATUS_CODE)
    IF (STATUS_CODE /= 0) THEN
        CALL log_error(ERR_MAIN_EXT_CMD)
        STOP ERR_MAIN_EXT_CMD
    END IF

    ! Continue with other code if STATUS_CODE is 0
    WRITE(*,*) 'PROCEEDING TO NEXT STEP.'

    ! --- Read Flowfield ---
    ! Read flowfield data
    CALL read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, &
                        Xgrid, Ygrid, Zgrid, data_count, error_status)
    
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_READ_FLOW)
        CALL cleanup_allocations()
        STOP ERR_MAIN_READ_FLOW
    END IF

    ! Generate initial disturbance
    CALL initial_disturbance(data_count, dist_mag, pert_0, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_DISTURBANCE)
        CALL cleanup_allocations()
        STOP ERR_MAIN_DISTURBANCE
    END IF

    ! Write flowfield data to file
    CALL write_flowfield_data(Xgrid, Ygrid, Zgrid, rho_in, p_in, T_in, &
                              U_in, V_in, W_in, data_count, error_status)
    IF (error_status /= 0) THEN
        CALL log_error(ERR_MAIN_WRITE_OUTPUT)
        CALL cleanup_allocations()
        STOP ERR_MAIN_WRITE_OUTPUT
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
        IF (ALLOCATED(pert_0)) DEALLOCATE(pert_0)
    END SUBROUTINE cleanup_allocations

END PROGRAM main
