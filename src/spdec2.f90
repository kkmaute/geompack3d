subroutine spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc,  &
  maxhv, maxpv, maxiw, maxwk, holv, vcl, regnum, hvl, pvl, iang, iwk, &
  wk, ierror )

!*****************************************************************************80
!
!! SPDEC2 decomposes a polygonal region with holes into simple polygons.
!
!  Discussion:
!
!    This routine decomposes a general polygonal region with interfaces and
!    holes into simple polygons using the vertex coordinate list,
!    head vertex list, and polygon vertex list data structures.
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
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of a
!    separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in radians
!    used in accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions
!    used in PVL array.
!
!    Input, integer ( kind = 4 ) NHOLE, the number of holes and hole interfaces.
!
!    Input, integer ( kind = 4 ) NHOLA, the number of 'attached' holes; these holes are
!    attached to the outer boundary of a subregion through vertices
!    or cut interfaces and have their edges in consecutive order on the
!    boundary.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array, should
!    be greater than or equal to the number of vertex coordinates required
!    for decomposition.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays,
!    should be greater than or equal to the number of polygons required
!    for decomposition.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays;
!    should be greater than or equal to the number of polygon vertices
!    required for decomposition.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    about 3 times maximum number of vertices in any polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    about 5 times maximum number of vertices in any polygon.
!
!    Input, integer ( kind = 4 ) HOLV(1:NHOLE*2+NHOLA), the indices in PVL of bottom or top
!    vertex of holes; first (next) NHOLE entries are for top (bottom)
!    vertices of holes and hole interfaces, with top (bottom)
!    vertices sorted in decreasing (increasing) lexicographic
!    (y,x) order of coord; last NHOLA entries are for attached
!    holes; if bottom vertex of attached hole is a simple
!    vertex of boundary curve containing the hole then entry
!    contains index of bottom vertex otherwise entry contains
!    index of top vertex (which is simple).
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles; see routine DSPGDC for more
!    details.  Note: The data structures should be as output from routines
!    DSMCPR or DSPGDC.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) angspc
  real    ( kind = 8 ) angtol
  logical              ci
  logical              cj
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) holv(*)
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) piptol
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real    ( kind = 8 ) wk(maxwk)

  tol = 100.0D+00 * epsilon ( tol )
!
!  For each simple hole, find cut edge from top vertex of hole to
!  a point on the outer boundary above top vertex, and update
!  VCL, HVL, PVL, IANG.
!
  piptol = pi + tol

  do i = 1, nhole

    call jnhole ( holv(i), angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
      maxwk, vcl, hvl, pvl, iang, iwk, wk, ierror )

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SPDEC2 - Fatal error!'
      write ( *, * ) '  JNHOLE returned an error condition.'
      return
    end if

  end do
!
!  Resolve remaining vertices in HOLV array if they are reflex
!  vertices. These vertices may no longer be reflex if they are the
!  endpoint of a cut edge from the top vertex of another hole or
!  of a previous separator.
!
  do i = nhole+1, nhole+nhole+nhola

    vr = holv(i)

    if ( piptol < iang(vr) ) then

      call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
        maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call insed2 ( vr, w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
        pvl, iang, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( 0 < w2 ) then

        call insed2 ( vr, w2, npolg, nvert, maxhv, maxpv, &
          vcl, regnum, hvl, pvl, iang, ierror )

        if ( ierror /= 0 ) then
          return
        end if

      end if

    end if

  end do

  if ( nhola == 0 ) then
    return
  end if
!
!  Check that polygons are simple. If polygon is simply-connected and
!  not simple then find a simple reflex vertex in polygon to resolve.
!
  p = 1

30 continue

  if ( npolg < p ) then
    return
  end if

  i = hvl(p)

  do

    if ( pvl(polg,pvl(edgv,i)) == p ) then
      go to 50
    end if

    i = pvl(succ,i)

    if ( i == hvl(p) ) then
      exit
    end if

  end do

  p = p + 1
  go to 30

50 continue

  ci = .true.

  do

    j = pvl(succ,i)
    cj = ( pvl(polg,pvl(edgv,j)) == p )

    if ( .not. ci .and. .not. cj .and. piptol < iang(j) ) then
      exit
    end if

    i = j
    ci = cj

  end do

  vr = j
  call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
    maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call insed2 ( vr, w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
    pvl, iang, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  if ( 0 < w2 ) then

    call insed2 ( vr, w2, npolg, nvert, maxhv, maxpv, &
      vcl, regnum, hvl, pvl, iang, ierror )

    if ( ierror /= 0 ) then
      return
    end if

  end if

  go to 30

end
