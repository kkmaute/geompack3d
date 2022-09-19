subroutine frstet ( shift, nv, vcl, map, i3, i4, ierr )

!*****************************************************************************80
!
!! FRSTET shifts vertices so the first 4 vertices are in general position in 3D.
!
!  Discussion:
!
!    This routine shifts or swaps vertices if necessary so first 3 vertices
!    (according to MAP) are not collinear and first 4 vertices are
!    not coplanar (so that first tetrahedron is valid).
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
!    Input, logical SHIFT, if TRUE, MAP(3), MAP(4) may be updated due to shift,
!    else they may be updated due to swaps; in former case,
!    it is assumed MAP gives vertices in lexicographic order.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) MAP(1:NV), on input, contains vertex indices of VCL.
!    On output, shifted or 2 swaps applied if necessary so that vertices
!    indexed by MAP(1), MAP(2), MAP(3), MAP(4) not coplanar.
!
!    Output, integer ( kind = 4 ) I3, I4, the indices such that MAP_in(I3) = MAP_out(3) and
!    MAP_in(I4) = MAP_out(4).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nv

  real    ( kind = 8 ) cmax
  real    ( kind = 8 ) cp1
  real    ( kind = 8 ) cp2
  real    ( kind = 8 ) cp3
  real    ( kind = 8 ) dmax
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dv2(3)
  real    ( kind = 8 ) dvk(3)
  real    ( kind = 8 ) dvl(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) map(nv)
  logical              shift
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,*)
!
!  First check that consecutive vertices are not identical.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( shift ) then
    l = nv - 1
  else
    l = 1
  end if

  m1 = map(1)

  do i = 1, l

    m = m1
    m1 = map(i+1)

    do k = 1, 3
      cmax = max ( abs ( vcl(k,m) ), abs ( vcl(k,m1) ) )
      if ( tol * cmax < abs ( vcl(k,m) - vcl(k,m1) ) .and. tol < cmax ) then
        go to 20
      end if
    end do

    ierr = 302
    return

20  continue

  end do
!
!  Find index K = I3 and L = I4.
!
  m1 = map(1)
  m2 = map(2)

  dv2(1:3) = vcl(1:3,m2) - vcl(1:3,m1)

  cmax = max ( abs ( vcl(1,m1) ), abs ( vcl(2,m1) ), abs ( vcl(3,m1) ), &
    abs ( vcl(1,m2) ), abs ( vcl(2,m2) ), abs ( vcl(3,m2) ) )
  k = 2

  do

    k = k + 1

    if ( nv < k ) then
      ierr = 303
      return
    end if

    m = map(k)

    dvk(1:3) = vcl(1:3,m) - vcl(1:3,m1)

    dmax = max ( cmax, abs ( vcl(1,m) ), abs ( vcl(2,m) ), abs ( vcl(3,m) ) )

    cp1 = dv2(2) * dvk(3) - dv2(3) * dvk(2)
    cp2 = dv2(3) * dvk(1) - dv2(1) * dvk(3)
    cp3 = dv2(1) * dvk(2) - dv2(2) * dvk(1)

    if ( tol * dmax < max ( abs ( cp1 ), abs ( cp2 ), abs ( cp3 ) ) ) then
      exit
    end if

  end do

  cmax = dmax
  l = k

  do

    l = l + 1

    if ( nv < l ) then
      ierr = 304
      return
    end if

    m = map(l)

    dvl(1:3) = vcl(1:3,m) - vcl(1:3,m1)

    dmax = max ( cmax, abs ( vcl(1,m) ), abs ( vcl(2,m) ), abs ( vcl(3,m) ) )

    dotp = dvl(1) * cp1 + dvl(2) * cp2 + dvl(3) * cp3

    if ( tol * dmax < abs ( dotp ) ) then
      exit
    end if

  end do
!
!  Shift or swap elements of MAP if necessary.
!
  if ( shift ) then

    if ( 3 < k ) then
      m1 = map(k)
    end if

    if ( 4 < l ) then
      m2 = map(l)
      do i = l, k+2, -1
         map(i) = map(i-1)
      end do
      do i = k+1, 5, -1
        map(i) = map(i-2)
      end do
      map(4) = m2
    end if

    if ( 3 < k ) then
      map(3) = m1
    end if

  else

    if ( 3 < k ) then
      m = map(3)
      map(3) = map(k)
      map(k) = m
    end if

    if ( 4 < l ) then
      m = map(4)
      map(4) = map(l)
      map(l) = m
    end if

  end if

  i3 = k
  i4 = l

  return
end
