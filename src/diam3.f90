subroutine diam3 ( nvrt, vcl, i1, i2, diamsq )

!*****************************************************************************80
!
!! DIAM3 finds the diameter of a set of 3D points.
!
!  Discussion:
!
!    This routine computes the diameter (largest distance) of set of 3D points,
!    and returns two vertex indices realizing diameter.
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
!    Input, integer ( kind = 4 ) NVRT, the number of vertices.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NVRT), the vertex coordinate list.
!
!    Output, integer ( kind = 4 ) I1, I2, the vertex indices realizing diameter (I1 < I2).
!
!    Output, real ( kind = 8 ) DIAMSQ, the square of diameter.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) diamsq
  real    ( kind = 8 ) distsq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  real    ( kind = 8 ) vcl(3,nvrt)

  diamsq = -1.0D+00

  do i = 1, nvrt-1

    do j = i+1, nvrt

      distsq = (vcl(1,i) - vcl(1,j))**2 + (vcl(2,i) - vcl(2,j))**2 &
        + (vcl(3,i) - vcl(3,j))**2

      if ( diamsq < distsq ) then
        diamsq = distsq
        i1 = i
        i2 = j
      end if

    end do

  end do

  return
end
