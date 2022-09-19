subroutine width2 ( nvrt, xc, yc, i1, i2, widsq, ierr )

!*****************************************************************************80
!
!! WIDTH2 determines the width of a convex polygon.
!
!  Discussion:
!
!    This routine finds the width (minimum breadth) of a convex polygon with
!    vertices given in counterclockwise order and with all interior
!    angles < PI.
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
!    Input, integter NVRT, number of vertices on the boundary of convex
!    polygon.
!
!    Input, real ( kind = 8 ) XC(1:NVRT), YC(1:NVRT), the vertex coordinates
!    in counterclockwise order.
!
!    Output, integer ( kind = 4 ) I1, I2, the indices in XC,YC such that width is from vertex
!    (XC(I1),YC(I1)) to line joining (XC(I2),YC(I2)) and
!    (XC(I2+1),YC(I2+1)), where index NVRT+1 is same as 1.
!
!    Output, real ( kind = 8 ) WIDSQ, the square of width.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) a
  real    ( kind = 8 ) area1
  real    ( kind = 8 ) area2
  real    ( kind = 8 ) areatr
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  real    ( kind = 8 ) c1mtol
  real    ( kind = 8 ) c1ptol
  real    ( kind = 8 ) dist
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  logical              first
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) m
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) widsq
  real    ( kind = 8 ) xc(nvrt)
  real    ( kind = 8 ) yc(nvrt)
!
!  Find first vertex which is farthest from edge connecting
!  vertices with indices NVRT, 1.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  first = .true.
  c1mtol = 1.0D+00 - tol
  c1ptol = 1.0D+00 + tol
  j = nvrt
  jp1 = 1
  k = 2
  area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )

10 continue

  area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

  if ( area1 * c1ptol < area2 ) then
    area1 = area2
    k = k + 1
    go to 10
  end if

  m = k
  widsq = 0.0D+00
!
!  Find width = min distance of antipodal edge-vertex pairs.
!
20 continue

  kp1 = k + 1
  if ( nvrt < kp1 ) then
    kp1 = 1
  end if

  area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(kp1), yc(kp1) )

  if ( area1 * c1ptol < area2 ) then
    a = j
    b = k
    k = k + 1
    c = k
    if ( nvrt < c ) then
      c = 1
    end if
    area1 = area2
  else if ( area2 < area1*c1mtol ) then
    a = k
    b = j
    c = jp1
    j = jp1
    jp1 = j + 1
    area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
  else
    a = k
    b = j
    c = jp1
    k = k + 1
    j = jp1
    jp1 = j + 1
    area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
  end if

  if ( m < j .or. nvrt < k ) then

    if ( first .and. 2 < m ) then
!
!  Possibly restart with M decreased by 1.
!
      first = .false.
      m = m - 1
      j = nvrt
      jp1 = 1
      k = m
      area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
      area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

      if ( area2 <= area1 * ( 1.0D+00 + 25.0D+00 * tol ) ) then
        area1 = area2
        widsq = 0.0D+00
        go to 20
      end if

    end if

    ierr = 201
    return

  end if

  dx = xc(c) - xc(b)
  dy = yc(c) - yc(b)
  dist = ((yc(a) - yc(b))*dx - (xc(a) - xc(b))*dy)**2/ (dx**2 + dy**2)

  if ( dist < widsq .or. widsq <= 0.0D+00 ) then
    widsq = dist
    i1 = a
    i2 = b
  end if

  if ( j /= m .or. k /= nvrt) go to 20

  return
end
