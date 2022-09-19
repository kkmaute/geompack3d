function umdf3 ( x, y, z )

!*****************************************************************************80
!
!! UMDF3 is a sample user mesh distribution function for 3D.
!
!  Discussion:
!
!    This routine is an example of a user-supplied mesh distribution function
!    which may be used in place of a heuristic mesh distribution function, or
!    may be a dummy routine.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, Z, the coordinates of 3D point
!
!    Output, real ( kind = 8 ) UMDF3, the mesh distribution function
!    value at (X,Y,Z)
!
  implicit none

  real    ( kind = 8 ) umdf3
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z
!
!  Mesh distribution function for airplane exterior region
!
  if ( -1.5D+00 <= x .and. x <= 38.0D+00 .and. &
                      abs(y) <= 4.0D+00  .and. &
                      abs(z) <= 4.0D+00 ) then

     umdf3 = 8.0D+00

  else if ( 11.0D+00 <= x .and. x <= 18.0D+00 .and. &
            -2.8D+00 <= y .and. y <=  2.0D+00 .and. &
                           abs(z) <= 16.5D+00 ) then

     umdf3 = 8.0D+00

  else if ( 31.0D+00 <= x .and. x <= 37.0D+00 .and. &
            -2.5D+00 <= y .and. y <=  2.0D+00 .and. &
                           abs(z) <= 9.0D+00 ) then

     umdf3 = 8.0D+00

  else if ( 31.0D+00 <= x .and. x <= 37.0D+00 .and. &
             0.0D+00 <= y .and. y <=  9.5D+00 .and. &
                           abs(z) <= 4.0D+00 ) then

     umdf3 = 8.0D+00

  else if ( -3.0D+00 <= x .and. x <= 40.0D+00 .and. &
            -9.0D+00 <= y .and. y <= 11.0D+00 .and. &
                           abs(z) <= 18.0D+00 ) then
     umdf3 = 3.0D+00
  else
     umdf3 = 1.0D+00
  end if

  return
end
