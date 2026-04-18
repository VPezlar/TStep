! =============================================================================
! error_handling  --  centralized error codes + uniform logging.
!
! Every numbered error in TStep is defined here as an INTEGER(ik) PARAMETER,
! grouped by originating module in blocks of 100:
!   000-099 Main          400-499 Flow reading
!   100-199 Setup         500-599 Output writing
!   200-299 OpenFOAM scalar I/O   600-699 External commands
!   300-399 OpenFOAM vector I/O   700-799 Random disturbance
!                                 800-899 Arnoldi
!                                 900-999 Eigendata writing
!
! `log_error(code [, info])` prints a banner with the code, module name, and
! human-readable description to stdout; `get_module_name` and
! `get_error_description` provide the text. Adding a new error means: add the
! PARAMETER, export it from the PRIVATE/PUBLIC list below, and add a CASE in
! `get_error_description`.
! =============================================================================
MODULE error_handling
    USE accuracy
    
    IMPLICIT NONE
    
    ! ============================================================
    ! CENTRALIZED ERROR CODE DEFINITIONS
    ! ============================================================
    ! Structure: ERROR_CODE = MODULE_ID * 100 + SPECIFIC_ERROR
    ! Each module owns a range of 100 error codes
    
    ! --- Main Program (0-99) ---
    INTEGER(ik), PARAMETER :: ERR_MAIN_CONFIG = 1
    INTEGER(ik), PARAMETER :: ERR_MAIN_EXT_CMD = 2
    INTEGER(ik), PARAMETER :: ERR_MAIN_READ_FLOW = 3
    INTEGER(ik), PARAMETER :: ERR_MAIN_DISTURBANCE = 4
    INTEGER(ik), PARAMETER :: ERR_MAIN_WRITE_OUTPUT = 5
    
    ! --- Setup Module (100-199) ---
    INTEGER(ik), PARAMETER :: ERR_SETUP_FILE_OPEN = 101
    INTEGER(ik), PARAMETER :: ERR_SETUP_NAMELIST_READ = 102
    INTEGER(ik), PARAMETER :: ERR_SETUP_INVALID_PARAM = 103
    
    ! --- OpenFOAM I/O - Scalars (200-299) ---
    INTEGER(ik), PARAMETER :: ERR_SCALAR_OPEN = 201
    INTEGER(ik), PARAMETER :: ERR_SCALAR_HEADER = 202
    INTEGER(ik), PARAMETER :: ERR_SCALAR_COUNT = 203
    INTEGER(ik), PARAMETER :: ERR_SCALAR_INVALID_COUNT = 204
    INTEGER(ik), PARAMETER :: ERR_SCALAR_ALLOC = 205
    INTEGER(ik), PARAMETER :: ERR_SCALAR_SKIP = 206
    INTEGER(ik), PARAMETER :: ERR_SCALAR_READ = 207
    
    ! --- OpenFOAM I/O - Vectors (300-399) ---
    INTEGER(ik), PARAMETER :: ERR_VECTOR_OPEN = 301
    INTEGER(ik), PARAMETER :: ERR_VECTOR_HEADER = 302
    INTEGER(ik), PARAMETER :: ERR_VECTOR_COUNT = 303
    INTEGER(ik), PARAMETER :: ERR_VECTOR_SKIP = 304
    INTEGER(ik), PARAMETER :: ERR_VECTOR_INVALID_COUNT = 305
    INTEGER(ik), PARAMETER :: ERR_VECTOR_ALLOC = 306
    INTEGER(ik), PARAMETER :: ERR_VECTOR_EOF = 307
    INTEGER(ik), PARAMETER :: ERR_VECTOR_FORMAT = 308
    INTEGER(ik), PARAMETER :: ERR_VECTOR_READ = 309
    
    ! --- Flow Reading Module (400-499) ---
    INTEGER(ik), PARAMETER :: ERR_FLOW_PRESSURE = 401
    INTEGER(ik), PARAMETER :: ERR_FLOW_DENSITY = 402
    INTEGER(ik), PARAMETER :: ERR_FLOW_TEMPERATURE = 403
    INTEGER(ik), PARAMETER :: ERR_FLOW_VELOCITY = 404
    INTEGER(ik), PARAMETER :: ERR_FLOW_GRID = 405
    INTEGER(ik), PARAMETER :: ERR_FLOW_UNKNOWN_FORMAT = 406
    
    ! --- Output Writing Module (500-599) ---
    INTEGER(ik), PARAMETER :: ERR_OUTPUT_FILE_OPEN = 501
    INTEGER(ik), PARAMETER :: ERR_OUTPUT_WRITE_HEADER = 502
    INTEGER(ik), PARAMETER :: ERR_OUTPUT_WRITE_DATA = 503
    
    ! --- External Commands Module (600-699) ---
    INTEGER(ik), PARAMETER :: ERR_CMD_LAUNCH_FAILED = 601
    INTEGER(ik), PARAMETER :: ERR_CMD_NONZERO_EXIT = 602
    
    ! --- Random Disturbance Module (700-799) ---
    INTEGER(ik), PARAMETER :: ERR_DIST_INVALID_LENGTH = 701
    INTEGER(ik), PARAMETER :: ERR_DIST_ALLOC = 702
    INTEGER(ik), PARAMETER :: ERR_DIST_ZERO_NORM = 703
    
    ! --- Arnoldi Module (800-899) ---
    INTEGER(ik), PARAMETER :: ERR_ARNOLDI_ZERO_V1 = 801
    INTEGER(ik), PARAMETER :: ERR_ARNOLDI_INVALID_DIM = 802
    INTEGER(ik), PARAMETER :: ERR_ARNOLDI_INVALID_KRYLOV = 803
    INTEGER(ik), PARAMETER :: ERR_ARNOLDI_ALLOC = 804
    INTEGER(ik), PARAMETER :: ERR_ARNOLDI_LAPACK = 805
    INTEGER(ik), PARAMETER :: ERR_ARNOLDI_NOT_IMPLEMENTED = 806
    
    ! --- Eigendata Writing Module (900-999) ---
    INTEGER(ik), PARAMETER :: ERR_EIGENDATA_OPEN_EVAL = 901
    INTEGER(ik), PARAMETER :: ERR_EIGENDATA_OPEN_EVEC = 902
    INTEGER(ik), PARAMETER :: ERR_EIGENDATA_WRITE_EVAL = 903
    INTEGER(ik), PARAMETER :: ERR_EIGENDATA_WRITE_EVEC = 904
    
    ! ============================================================
    ! ERROR MESSAGE LOOKUP
    ! ============================================================
    
    PRIVATE
    PUBLIC :: ERR_MAIN_CONFIG, ERR_MAIN_EXT_CMD, ERR_MAIN_READ_FLOW, &
              ERR_MAIN_DISTURBANCE, ERR_MAIN_WRITE_OUTPUT, &
              ERR_SETUP_FILE_OPEN, ERR_SETUP_NAMELIST_READ, ERR_SETUP_INVALID_PARAM, &
              ERR_SCALAR_OPEN, ERR_SCALAR_HEADER, ERR_SCALAR_COUNT, &
              ERR_SCALAR_INVALID_COUNT, ERR_SCALAR_ALLOC, ERR_SCALAR_SKIP, ERR_SCALAR_READ, &
              ERR_VECTOR_OPEN, ERR_VECTOR_HEADER, ERR_VECTOR_COUNT, ERR_VECTOR_SKIP, &
              ERR_VECTOR_INVALID_COUNT, ERR_VECTOR_ALLOC, ERR_VECTOR_EOF, ERR_VECTOR_FORMAT, &
              ERR_VECTOR_READ, &
              ERR_FLOW_PRESSURE, ERR_FLOW_DENSITY, ERR_FLOW_TEMPERATURE, &
              ERR_FLOW_VELOCITY, ERR_FLOW_GRID, ERR_FLOW_UNKNOWN_FORMAT, &
              ERR_OUTPUT_FILE_OPEN, ERR_OUTPUT_WRITE_HEADER, ERR_OUTPUT_WRITE_DATA, &
              ERR_CMD_LAUNCH_FAILED, ERR_CMD_NONZERO_EXIT, &
              ERR_DIST_INVALID_LENGTH, ERR_DIST_ALLOC, ERR_DIST_ZERO_NORM, &
              ERR_ARNOLDI_ZERO_V1, ERR_ARNOLDI_INVALID_DIM, ERR_ARNOLDI_INVALID_KRYLOV, &
              ERR_ARNOLDI_ALLOC, ERR_ARNOLDI_LAPACK, ERR_ARNOLDI_NOT_IMPLEMENTED, &
              ERR_EIGENDATA_OPEN_EVAL, ERR_EIGENDATA_OPEN_EVEC, &
              ERR_EIGENDATA_WRITE_EVAL, ERR_EIGENDATA_WRITE_EVEC, &
              log_error, get_module_name, get_error_description
    
CONTAINS

    ! Get module name from error code
    FUNCTION get_module_name(error_code) RESULT(module_name)
        INTEGER(ik), INTENT(in) :: error_code
        CHARACTER(len=30) :: module_name
        INTEGER(ik) :: module_id
        
        module_id = error_code / 100
        
        SELECT CASE (module_id)
            CASE (0)
                module_name = "MAIN"
            CASE (1)
                module_name = "SETUP"
            CASE (2)
                module_name = "OPENFOAM_IO_SCALAR"
            CASE (3)
                module_name = "OPENFOAM_IO_VECTOR"
            CASE (4)
                module_name = "FLOW_READING"
            CASE (5)
                module_name = "OUTPUT_WRITING"
            CASE (6)
                module_name = "EXTERNAL_COMMAND"
            CASE (7)
                module_name = "RANDOM_DISTURBANCE"
            CASE (8)
                module_name = "ARNOLDI"
            CASE (9)
                module_name = "EIGENDATA_WRITING"
            CASE DEFAULT
                module_name = "UNKNOWN_MODULE"
        END SELECT
        
    END FUNCTION get_module_name
    
    ! Get human-readable error description
    FUNCTION get_error_description(error_code) RESULT(description)
        INTEGER(ik), INTENT(in) :: error_code
        CHARACTER(len=100) :: description
        
        SELECT CASE (error_code)
            ! Main program errors
            CASE (ERR_MAIN_CONFIG)
                description = "Configuration read failed"
            CASE (ERR_MAIN_EXT_CMD)
                description = "External command execution failed"
            CASE (ERR_MAIN_READ_FLOW)
                description = "Flow field read failed"
            CASE (ERR_MAIN_DISTURBANCE)
                description = "Disturbance generation failed"
            CASE (ERR_MAIN_WRITE_OUTPUT)
                description = "Output write failed"
                
            ! Setup errors
            CASE (ERR_SETUP_FILE_OPEN)
                description = "Cannot open configuration file"
            CASE (ERR_SETUP_NAMELIST_READ)
                description = "Cannot read namelist"
            CASE (ERR_SETUP_INVALID_PARAM)
                description = "Invalid parameter value"
                
            ! Scalar I/O errors
            CASE (ERR_SCALAR_OPEN)
                description = "Cannot open scalar file"
            CASE (ERR_SCALAR_HEADER)
                description = "Error reading scalar file header"
            CASE (ERR_SCALAR_COUNT)
                description = "Cannot read scalar data count"
            CASE (ERR_SCALAR_INVALID_COUNT)
                description = "Invalid scalar data count"
            CASE (ERR_SCALAR_ALLOC)
                description = "Memory allocation failed for scalars"
            CASE (ERR_SCALAR_SKIP)
                description = "Error skipping scalar discard line"
            CASE (ERR_SCALAR_READ)
                description = "Error reading scalar data"
                
            ! Vector I/O errors
            CASE (ERR_VECTOR_OPEN)
                description = "Cannot open vector file"
            CASE (ERR_VECTOR_HEADER)
                description = "Error reading vector file header"
            CASE (ERR_VECTOR_COUNT)
                description = "Cannot read vector data count"
            CASE (ERR_VECTOR_SKIP)
                description = "Error skipping vector discard line"
            CASE (ERR_VECTOR_INVALID_COUNT)
                description = "Invalid vector data count"
            CASE (ERR_VECTOR_ALLOC)
                description = "Memory allocation failed for vectors"
            CASE (ERR_VECTOR_EOF)
                description = "Premature end of file in vector data"
            CASE (ERR_VECTOR_FORMAT)
                description = "Vector data format error"
            CASE (ERR_VECTOR_READ)
                description = "Error reading vector data"
                
            ! Flow reading errors
            CASE (ERR_FLOW_PRESSURE)
                description = "Failed to read pressure field"
            CASE (ERR_FLOW_DENSITY)
                description = "Failed to read density field"
            CASE (ERR_FLOW_TEMPERATURE)
                description = "Failed to read temperature field"
            CASE (ERR_FLOW_VELOCITY)
                description = "Failed to read velocity field"
            CASE (ERR_FLOW_GRID)
                description = "Failed to read grid coordinates"
            CASE (ERR_FLOW_UNKNOWN_FORMAT)
                description = "Unknown flow format specified"
                
            ! Output writing errors
            CASE (ERR_OUTPUT_FILE_OPEN)
                description = "Cannot open output file"
            CASE (ERR_OUTPUT_WRITE_HEADER)
                description = "Error writing output header"
            CASE (ERR_OUTPUT_WRITE_DATA)
                description = "Error writing output data"
                
            ! External command errors
            CASE (ERR_CMD_LAUNCH_FAILED)
                description = "Command failed to launch"
            CASE (ERR_CMD_NONZERO_EXIT)
                description = "Command returned non-zero exit code"
                
            ! Random disturbance errors
            CASE (ERR_DIST_INVALID_LENGTH)
                description = "Invalid disturbance vector length"
            CASE (ERR_DIST_ALLOC)
                description = "Memory allocation failed for disturbance"
            CASE (ERR_DIST_ZERO_NORM)
                description = "Disturbance vector has zero norm"
                
            ! Arnoldi errors
            CASE (ERR_ARNOLDI_ZERO_V1)
                description = "Initial Arnoldi vector has zero norm"
            CASE (ERR_ARNOLDI_INVALID_DIM)
                description = "Invalid matrix dimension for Arnoldi"
            CASE (ERR_ARNOLDI_INVALID_KRYLOV)
                description = "Invalid Krylov subspace size"
            CASE (ERR_ARNOLDI_ALLOC)
                description = "Memory allocation failed in Arnoldi"
            CASE (ERR_ARNOLDI_LAPACK)
                description = "LAPACK eigenvalue computation failed"
            CASE (ERR_ARNOLDI_NOT_IMPLEMENTED)
                description = "Arnoldi matvec (Frechet derivative) not yet implemented"
                
            ! Eigendata writing errors
            CASE (ERR_EIGENDATA_OPEN_EVAL)
                description = "Cannot open eigenvalues.dat for writing"
            CASE (ERR_EIGENDATA_OPEN_EVEC)
                description = "Cannot open eigenvectors.dat for writing"
            CASE (ERR_EIGENDATA_WRITE_EVAL)
                description = "Error writing eigenvalues data"
            CASE (ERR_EIGENDATA_WRITE_EVEC)
                description = "Error writing eigenvectors data"
                
            CASE DEFAULT
                description = "Unknown error code"
        END SELECT
        
    END FUNCTION get_error_description
    
    ! Unified error logging subroutine
    SUBROUTINE log_error(error_code, additional_info)
        INTEGER(ik), INTENT(in) :: error_code
        CHARACTER(len=*), INTENT(in), OPTIONAL :: additional_info
        
        CHARACTER(len=30) :: module_name
        CHARACTER(len=100) :: description
        
        module_name = get_module_name(error_code)
        description = get_error_description(error_code)
        
        WRITE(*,'(A)') '========================================='
        WRITE(*,'(A,I0)') 'ERROR CODE: ', error_code
        WRITE(*,'(A,A)') 'MODULE:     ', TRIM(module_name)
        WRITE(*,'(A,A)') 'DESCRIPTION: ', TRIM(description)
        IF (PRESENT(additional_info)) THEN
            WRITE(*,'(A,A)') 'DETAILS:    ', TRIM(additional_info)
        END IF
        WRITE(*,'(A)') '========================================='
        
    END SUBROUTINE log_error

END MODULE error_handling
