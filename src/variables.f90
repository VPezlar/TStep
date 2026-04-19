! =============================================================================
! variables  --  shared configuration globals populated from inputs/inputs.in.
!
! This is the central bag of runtime settings that every other module reads
! after `configurationRead` in setup.f90 has parsed the namelists. Keeping
! them here (rather than passing dozens of arguments) lets the solver-agnostic
! pipeline switch formats or parameters without plumbing changes.
!
! Groups:
!   - I/O paths + formats (file_grid_in, file_var_in, file_var_out, ...)
!   - Simulation driver (COMMAND_RUN, flow_format, N_HEADER_grid/var, dist_mag)
!   - Arnoldi parameters (krylov_size, frechet_order, eps_0, TTime, ...)
!   - Performance (num_threads)
! =============================================================================
MODULE variables

    USE accuracy

    IMPLICIT NONE

    INTEGER(ik) :: N_HEADER_grid, N_HEADER_var
    CHARACTER(len=256) :: file_grid_in, file_var_in, file_var_out, output_file, flow_format, COMMAND_RUN
    REAL(rk) :: dist_mag
    
    ! Arnoldi parameters
    INTEGER(ik) :: krylov_size
    INTEGER(ik) :: frechet_order
    CHARACTER(len=20) :: eigenvalue_sort_by
    REAL(rk) :: eps_0, TTime
    
    ! Performance parameters
    INTEGER(ik) :: num_threads  ! Number of threads for BLAS/LAPACK (0 = auto)

END MODULE variables
