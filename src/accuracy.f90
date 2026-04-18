! =============================================================================
! accuracy  --  project-wide numerical precision kinds.
!
! Defines two aliases every other module uses:
!   ik = int32   (default integer kind)
!   rk = real64  (default real/complex kind, i.e. double precision)
! Change these in one place here if the project ever needs a different
! precision; everything downstream follows automatically.
! =============================================================================
MODULE accuracy

    USE, INTRINSIC :: iso_fortran_env, ik=>int32, rk=>real64

END MODULE accuracy
