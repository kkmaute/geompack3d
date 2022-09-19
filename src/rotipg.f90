subroutine rotipg ( xeye, yeye, nvrt, xc, yc, ierr )

!*****************************************************************************80
!
!! ROTIPG rotates the indices of the vertices of a simple polygon.
!
!  Discussion:
!
!    This routine rotates the indices of the vertices of a simple polygon
!    and possibly inserts one vertex so that (XC(0),YC(0)) is the
!    point on the horizontal line through (XEYE,YEYE) and on the
!    boundary of the polygon which is closest to and to the right
!    of (XEYE,YEYE).  (XEYE,YEYE) is an eyepoint in the interior or
!    blocked exterior of the polygon.  In the former (latter) case,
!    the vertices must be in counterclockwise (CW) order.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of eyepoint.
!    On output, vertices are in same orientation, but with indices rotated,
!    and possibly (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)) has been added.
!
!    Input/output, integer ( kind = 4 ) NVRT, the number of vertices on boundary of simple
!    polygon.  On output, increased by 1 from input iff the closest vertex
!    is a new vertex.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertices of polygon
!    in counterclockwise (or clockwise) order if eyepoint is interior (or
!    blocked exterior); (XC(0),YC(0)) = (XC(NVRT),YC(NVRT))
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) irgt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) r
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) xc(0:nvrt+1)
  real    ( kind = 8 ) xeye
  real    ( kind = 8 ) xint
  real    ( kind = 8 ) xrgt
  real    ( kind = 8 ) xt
  real    ( kind = 8 ) yc(0:nvrt+1)
  real    ( kind = 8 ) yeye
  real    ( kind = 8 ) yeyemt
  real    ( kind = 8 ) yeyept
  real    ( kind = 8 ) yt

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  dy = 0.0D+00
  do i = 0, nvrt-1
    dy = max ( dy, abs ( yc(i+1) - yc(i) ) )
  end do

  yeyemt = yeye - tol * dy
  yeyept = yeye + tol * dy
  n = nvrt + 1
  irgt = n
  xrgt = 0.0D+00
!
!  Determine closest point on boundary which is to the right of
!  (XEYE,YEYE) and on the horizontal line through (XEYE,YEYE).
!  The closest point must be on an edge which intersects the
!  horizontal line and has (XEYE,YEYE) to its left.
!
  do i = 0, nvrt-1

    if ( yeyept < yc(i) .or. yc(i+1) < yeyemt ) then
      cycle
    end if

    if ( yc(i) < yeyemt .and. yeyept < yc(i+1) ) then

      xint = (yeye-yc(i))*(xc(i+1)-xc(i))/(yc(i+1)-yc(i)) + xc(i)

      if ( xeye < xint ) then
        if ( xint < xrgt .or. irgt == n ) then
          irgt = -(i + 1)
          xrgt = xint
        end if
      end if

    else if ( yeyemt <= yc(i) .and. yeyept < yc(i+1) ) then

      if ( xeye < xc(i) ) then
        if ( xc(i) < xrgt .or. irgt == n ) then
          irgt = i
          xrgt = xc(i)
        end if
      end if

    else if ( yc(i) < yeyemt .and. yc(i+1) <= yeyept ) then

      if ( xeye < xc(i+1) ) then
        if ( xc(i+1) < xrgt .or. irgt == n ) then
          irgt = i + 1
          xrgt = xc(i+1)
        end if
      end if
    end if

  end do

  if ( irgt == n ) then
    ierr = 205
    return
  end if

  if ( irgt == 0 .or. irgt == nvrt ) then
    return
  end if

  if ( irgt < 0 ) then

    irgt = -irgt
    do i = nvrt, irgt, -1
      xc(i+1) = xc(i)
      yc(i+1) = yc(i)
    end do

    xc(irgt) = xrgt
    yc(irgt) = yeye
    nvrt = nvrt + 1

  end if
!
!  Rotate the indices of the vertices so that (XC(IRGT),YC(IRGT))
!  becomes (XC(0),YC(0)). Compute A = GCD(NVRT,IRGT).
!
  a = nvrt
  b = irgt

  do

    r = mod ( a, b )
    a = b
    b = r

    if ( r <= 0 ) then
      exit
    end if

  end do

  m = nvrt / a - 1

  do i = 0, a-1

    xt = xc(i)
    yt = yc(i)
    k = i

    do j = 1, m
      l = k + irgt
      if ( nvrt <= l ) then
        l = l - nvrt
      end if
      xc(k) = xc(l)
      yc(k) = yc(l)
      k = l
    end do

    xc(k) = xt
    yc(k) = yt

  end do

  xc(nvrt) = xc(0)
  yc(nvrt) = yc(0)

  return
end
