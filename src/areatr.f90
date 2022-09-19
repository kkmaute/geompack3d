function areatr ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! AREATR computes twice the signed area of a triangle.
!
!  Discussion:
!
!    This routine computes twice the signed area of the triangle with
!    vertices (XA,YA), (XB,YB), and (XC,YC) in counterclockwise or
!    clockwise order.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XA, YA, XB, YB, XC, YC, the vertex coordinates.
!
!    Output, real ( kind = 8 ) AREATR, twice the signed area of triangle,
!    positive if counterclockwise.
!
  implicit none

  real    ( kind = 8 ) areatr
  real    ( kind = 8 ) xa
  real    ( kind = 8 ) xb
  real    ( kind = 8 ) xc
  real    ( kind = 8 ) ya
  real    ( kind = 8 ) yb
  real    ( kind = 8 ) yc

  areatr = ( xb - xa ) * ( yc - ya ) - ( xc - xa ) * ( yb - ya )

  return
end
