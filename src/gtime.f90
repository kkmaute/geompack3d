subroutine gtime ( time )

!*****************************************************************************80
!
!! GTIME returns the current CPU time in seconds.
!
!  Modified:
!
!    17 September 2001
!
!  Parameters:
!
!    Output, real TIME, the current CPU time in seconds.
!
  implicit none

  real    time

  call cpu_time ( time )

  return
end
