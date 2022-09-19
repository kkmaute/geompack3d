subroutine diam2 ( nvrt, xc, yc, i1, i2, diamsq, ierr )

!*****************************************************************************80
!
!! DIAM2 finds the diameter of a convex polygon.
!
!  Discussion:
!
!    This routine finds the diameter of a convex polygon with vertices
!    given in counterclockwise order and with all interior angles < PI.
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
!    convex polygon.
!
!    Input, real ( kind = 8 ) XC(1:NVRT), YC(1:NVRT), the vertex coordinates
!    in counterclockwise order.
!
!    Output, integer ( kind = 4 ) I1, I2, the indices in XC, YC of diameter edge; diameter
!    is from (XC(I1),YC(I1)) to (XC(I2),YC(I2)).
!
!    Output, real ( kind = 8 ) DIAMSQ, the square of diameter.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) area1
  real    ( kind = 8 ) area2
  real    ( kind = 8 ) areatr
  real    ( kind = 8 ) c1mtol
  real    ( kind = 8 ) c1ptol
  real    ( kind = 8 ) diamsq
  real    ( kind = 8 ) dist
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

  do

    area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k+1), yc(k+1) )

    if ( area2 <= area1 * c1ptol ) then
      exit
    end if

    area1 = area2
    k = k + 1

  end do

  m = k
  diamsq = 0.0D+00
!
!  Find diameter = maximum distance of antipodal pairs.
!
20 continue

  kp1 = k + 1

  if ( nvrt < kp1 ) then
    kp1 = 1
  end if

  area2 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(kp1), yc(kp1) )

  if ( area1 * c1ptol < area2 ) then
    k = k + 1
    area1 = area2
  else if ( area2 < area1 * c1mtol ) then
    j = jp1
    jp1 = j + 1
    area1 = areatr ( xc(j), yc(j), xc(jp1), yc(jp1), xc(k), yc(k) )
  else
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

      if ( area2 <= area1*(1.0D+00+25.0D+00 * tol) ) then
        area1 = area2
        diamsq = 0.0D+00
        go to 20
      end if

    end if

    ierr = 200
    return

  end if

  dist = ( xc(j) - xc(k) )**2 + ( yc(j) - yc(k) )**2

  if ( diamsq < dist ) then
    diamsq = dist
    i1 = j
    i2 = k
  end if

  if ( j /= m .or. k /= nvrt ) then
    go to 20
  end if

  return
end
