subroutine ptpolg ( dim, ldv, nv, inc, pgind, vcl, pt, nrml, dtol, inout )

!*****************************************************************************80
!
!! PTPOLG determines where a point lies with respect to a polygon.
!
!  Discussion:
!
!    This routine determines whether a point lies inside, outside, or on
!    boundary of a planar polygon in 2 or 3 dimensional space.
!    It is assumed that point lies in plane of polygon.
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
!    Input, integer ( kind = 4 ) DIM, the dimension of polygon (2 or 3).
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL array in calling routine.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices in polygon.
!
!    Input, integer ( kind = 4 ) INC, the increment for PGIND indicating indices of polygon.
!
!    Input, integer ( kind = 4 ) PGIND(0:NV*INC), the indices in VCL of polygon vertices
!    are in PGIND(0), PGIND(INC), ..., PGIND(NV*INC) with first and
!    last vertices identical.
!
!    Input, real ( kind = 8 ) VCL(1:DIM,1:*), the vertex coordinate list.
!
!    Input, real ( kind = 8 ) PT(1:DIM), the point for which in/out test
!    is applied.
!
!    Input, real ( kind = 8 ) NRML(1:3), the unit normal vector of plane
!    containing polygon, with vertices oriented counterclockwise with respect
!    to normal (used iff DIM = 3); normal is assumed to be (0,0,1) if DIM = 2.
!
!    Input, real ( kind = 8 ) DTOL, the absolute tolerance to determine
!    whether a point is on a line or plane.
!
!    Output, integer ( kind = 4 ) INOUT, +1, 0, or -1 depending on whether point PT is
!    inside polygon, on boundary of polygon, or outside polygon;
!    or -2 if error in input parameters.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) ldv

  real    ( kind = 8 ) area
  real    ( kind = 8 ) armax
  real    ( kind = 8 ) cp(3)
  real    ( kind = 8 ) da(3)
  real    ( kind = 8 ) db(3)
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) dist
  real    ( kind = 8 ) dtol
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) m
  real    ( kind = 8 ) nr(4)
  real    ( kind = 8 ) nrml(3)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) pgind(0:*)
  real    ( kind = 8 ) pt(dim)
  real    ( kind = 8 ) rhs(3)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) sa
  integer ( kind = 4 ) sb
  real    ( kind = 8 ) t
  real    ( kind = 8 ) ta
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(ldv,*)

  tol = 100.0D+00 * epsilon ( tol )
  inout = -2

  if ( dim < 2 .or. 3 < dim ) then
    return
  end if
!
!  Find edge subtending max area with PT as third triangle vertex.
!
  armax = 0.0D+00
  h = 0
  lb = pgind(0)
  db(1:dim) = vcl(1:dim,lb) - pt(1:dim)

  do i = 1, nv

    la = lb
    lb = pgind(i*inc)

    do j = 1, dim
      da(j) = db(j)
      db(j) = vcl(j,lb) - pt(j)
      dir(j) = vcl(j,lb) - vcl(j,la)
    end do

    if ( dim == 2 ) then
      area = abs ( da(1) * db(2) - db(1) * da(2) )
    else
      area = (da(2)*db(3)-db(2)*da(3))**2 + (da(3)*db(1)- &
        db(3)*da(1))**2 + (da(1)*db(2)-db(1)*da(2))**2
    end if

    if ( armax < area ) then
      h = i
      armax = area
    end if

  end do

  if ( dim == 2 ) then
    armax = armax**2
  end if

  if ( armax <= dtol**2 ) then
    return
  end if
!
!  Direction of ray is from PT through midpoint of edge subtending
!  max area. NR is unit normal of line or plane containing ray,
!  which is also orthogonal to NRML in 3D case.
!
  la = pgind((h-1)*inc)
  lb = pgind(h*inc)
  dir(1) = 0.5D+00*(vcl(1,la) + vcl(1,lb)) - pt(1)
  dir(2) = 0.5D+00*(vcl(2,la) + vcl(2,lb)) - pt(2)

  if ( dim == 2 ) then
    dist = sqrt(dir(1)**2 + dir(2)**2)
    dir(1) = dir(1) / dist
    dir(2) = dir(2) / dist
    dir(3) = 0.0D+00
    nr(1) = -dir(2)
    nr(2) = dir(1)
    nr(4) = nr(1)*pt(1) + nr(2)*pt(2)
    dist = nr(1)*vcl(1,lb) + nr(2)*vcl(2,lb) - nr(4)
  else
    dir(3) = 0.5D+00*(vcl(3,la) + vcl(3,lb)) - pt(3)
    dist = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)

    dir(1:3) = dir(1:3) / dist

    nr(1) = nrml(2)*dir(3) - nrml(3)*dir(2)
    nr(2) = nrml(3)*dir(1) - nrml(1)*dir(3)
    nr(3) = nrml(1)*dir(2) - nrml(2)*dir(1)

    nr(4) = dot_product ( nr(1:3), pt(1:3) )

    dist = dot_product ( nr(1:3), vcl(1:3,lb) ) - nr(4)

  end if

  if ( 0.0D+00 < dist ) then
    sb = 1
  else
    sb = -1
  end if

  m = 1
  if ( abs ( dir(1) ) < abs ( dir(2) ) ) then
    m = 2
  end if

  if ( abs ( dir(m) ) < abs ( dir(3) ) ) then
    m = 3
  end if

  k = 1
!
!  For remaining edges of polygon, check whether ray intersects edge.
!  Vertices or edges lying on ray are handled by looking at preceding
!  and succeeding edges not lying on ray.
!
  k = 1
  i = h + 1
  if ( nv < i ) then
    i = 1
  end if

40 continue

  la = lb
  lb = pgind(i*inc)
  sa = sb

  if ( dim == 2 ) then
    dist = nr(1)*vcl(1,lb) + nr(2)*vcl(2,lb) - nr(4)
  else
    dist = nr(1)*vcl(1,lb) + nr(2)*vcl(2,lb) + nr(3)*vcl(3,lb) - nr(4)
  end if

  if ( abs ( dist ) <= dtol ) then
    sb = 0
  else if ( 0.0D+00 < dist ) then
    sb = 1
  else
    sb = -1
  end if

  s = sa * sb

  if ( s < 0 ) then

    if ( dim == 2 ) then

      da(1) = vcl(1,la) - vcl(1,lb)
      da(2) = vcl(2,la) - vcl(2,lb)
      rhs(1) = vcl(1,la) - pt(1)
      rhs(2) = vcl(2,la) - pt(2)
      t = (rhs(1)*da(2) - rhs(2)*da(1))/(dir(1)*da(2) - dir(2)*da(1))

    else

      da(1:3) = vcl(1:3,la) - vcl(1:3,lb)
      rhs(1:3) = vcl(1:3,la) - pt(1:3)

      cp(1) = dir(2)*da(3) - dir(3)*da(2)
      cp(2) = dir(3)*da(1) - dir(1)*da(3)
      cp(3) = dir(1)*da(2) - dir(2)*da(1)

      l = 1
      if ( abs ( cp(1) ) < abs ( cp(2) ) ) then
        l = 2
      end if

      if ( abs ( cp(l) ) < abs ( cp(3) ) ) then
        l = 3
      end if

      if ( l == 1 ) then
        t = ( rhs(2) * da(3) - rhs(3) * da(2) ) / cp(1)
      else if ( l == 2 ) then
        t = ( rhs(3) * da(1) - rhs(1) * da(3) ) / cp(2)
      else
        t = ( rhs(1) * da(2) - rhs(2) * da(1) ) / cp(3)
      end if

    end if

    if ( dtol < t ) then
      k = k + 1
    else if ( -dtol <= t ) then
      inout = 0
      return
    end if

  else if ( s == 0 ) then

    l = lb

50  continue

    i = i + 1
    if ( nv < i ) then
      i = 1
    end if

    if ( i == h) then
      return
    end if

    la = lb
    lb = pgind(i*inc)

    if ( dim == 2 ) then
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
    else
      dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) + nr(3) * vcl(3,lb) - nr(4)
    end if

    if ( abs ( dist ) <= dtol ) then
      go to 50
    else if ( 0.0D+00 < dist ) then
      sb = 1
    else
      sb = -1
    end if

    t = ( vcl(m,l) - pt(m) ) / dir(m)

    if ( abs ( t ) <= dtol ) then
      inout = 0
      return
    end if

    if ( la /= l ) then
      ta = ( vcl(m,la) - pt(m) ) / dir(m)
      if ( abs ( ta ) <= dtol .or. t * ta < 0.0D+00 ) then
        inout = 0
        return
      end if
    end if

    if ( sa * sb < 0 .and. 0.0D+00 < t ) then
      k = k + 1
    end if

  end if

  i = i + 1
  if ( nv < i ) then
    i = 1
  end if

  if ( i /= h ) then
    go to 40
  end if
!
!  Point lies inside polygon if number of intersections K is odd.
!
  if ( mod ( k, 2 ) == 1 ) then
    inout = 1
  else
    inout = -1
  end if

  return
end
