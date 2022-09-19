function umdf2 ( x, y )

!*****************************************************************************80
!
!! UMDF2 is a sample user mesh distribution function for 2D.
!
!  Discussion:
!
!    This routine is a dummy user-supplied mesh distribution function which
!    is provided if a heuristic mesh distribution function is used.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of 2D point.
!
!    Output, real ( kind = 8 ) UMDF2, the mesh distribution function
!    value at (X,Y).
!
  implicit none

  real    ( kind = 8 ) umdf2
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  umdf2 = 1.0D+00

  return
end
