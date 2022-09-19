function areapg ( nvrt, xc, yc )

!*****************************************************************************80
!
!! AREAPG computes twice the signed area of a simple polygon.
!
!  Discussion:
!
!    This routine computes twice the signed area of a simple polygon,
!    with vertices given in circular (counterclockwise or clockwise) order.
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
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    polygon.  3 <= NVRT.
!
!    Input, real ( kind = 8 ) XC(1:NVRT), YC(1:NVRT), the vertex coordinates
!    in counterclockwise or clockwise order.
!
!    Output, real ( kind = 8 ) AREAPG, twice the signed area of the polygon,
!    positive if counterclockwise.
!
  implicit none

  real    ( kind = 8 ) areapg
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) i
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) xc(nvrt)
  real    ( kind = 8 ) yc(nvrt)

  sum2 = xc(1) * ( yc(2) - yc(nvrt) ) + xc(nvrt) * ( yc(1) - yc(nvrt-1) )
  do i = 2, nvrt-1
    sum2 = sum2 + xc(i) * ( yc(i+1) - yc(i-1) )
  end do

  areapg = sum2

  return
end
