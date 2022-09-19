subroutine rotpg ( nvrt, xc, yc, i1, i2, ibot, costh, sinth )

!*****************************************************************************80
!
!! ROTPG rotates a convex polygon.
!
!  Discussion:
!
!    This routine rotates a convex polygon so that a line segment joining two
!    of its vertices is parallel to the Y axis.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    the convex polygon.
!
!    Input/output, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT).  The vertex
!    coordinates in counterclockwise order;
!      (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!    On output, the rotated vertex coordinates; indices are
!    also rotated so that (XC(0),YC(0)) = (XC(NVRT),YC(NVRT))
!    is top vertex and (XC(IBOT),YC(IBOT)) is bottom vertex.
!
!    Input, integer ( kind = 4 ) I1, I2, the index of vertices of line segment;
!    I1, I2 must be greater than 0.
!
!    Output, integer ( kind = 4 ) IBOT, the index of bottom vertex.
!
!    Output, real ( kind = 8 ) COSTH, SINTH, the values COS(THETA) and
!    SIN(THETA) where THETA in [-PI,PI] is the rotation angle.
!
  implicit none

  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) costh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real    ( kind = 8 ) sinth
  real    ( kind = 8 ) theta
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) yc(0:nvrt)

  tol = 100.0D+00 * epsilon ( tol )

  itop = i1
  ibot = i2

  if ( yc(i1) == yc(i2) ) then
    if ( xc(i1) < xc(i2) ) then
      theta = - pi / 2.0D+00
    else
      theta = pi / 2.0D+00
    end if
  else
    if ( yc(i1) < yc(i2) ) then
      itop = i2
      ibot = i1
    end if
    theta = pi / 2.0D+00 &
      - atan2 ( yc(itop) - yc(ibot), xc(itop) - xc(ibot) )
  end if

  costh = cos(theta)
  sinth = sin(theta)

  do i = 1, nvrt
    x0 = xc(i)
    xc(i) = costh * x0 - sinth * yc(i)
    yc(i) = sinth * x0 + costh * yc(i)
  end do
!
!  Rotate indices.
!
  if ( itop /= nvrt ) then

    a = nvrt
    b = itop

    do

      r = mod ( a, b )
      a = b
      b = r

      if ( r <= 0 ) then
        exit
      end if

    end do

    m = nvrt / a - 1

    do i = 1, a

      x0 = xc(i)
      y0 = yc(i)
      k = i

      do j = 1, m
        l = k + itop
        if ( nvrt < l ) then
          l = l - nvrt
        end if
        xc(k) = xc(l)
        yc(k) = yc(l)
        k = l
      end do

      xc(k) = x0
      yc(k) = y0

    end do

    ibot = ibot - itop
    if ( ibot < 0 ) then
      ibot = ibot + nvrt
    end if

  end if

  xc(0) = xc(nvrt)
  yc(0) = yc(nvrt)

  return
end
