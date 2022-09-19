subroutine cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, &
  maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

!*****************************************************************************80
!
!! CVDEC2 decomposes a polygonal region into convex polygons.
!
!  Discussion:
!
!    This routine decomposes a general polygonal region (which is decomposed
!    into simple polygons on input) into convex polygons using
!    vertex coordinate list, head vertex list, and polygon vertex
!    list data structures.
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
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of
!    a separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in radians
!    used in accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or positions
!    used in PVL array.
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
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles; see routine DSPGDC for
!    more details.  Note that the data structures should be as output from
!    routine SPDEC2.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) angspc
  real    ( kind = 8 ) angtol
  integer ( kind = 4 ) hvl(maxhv)
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) piptol
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) v
  real    ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real    ( kind = 8 ) wk(maxwk)

  tol = 100.0D+00 * epsilon ( tol )
  ierror = 0
!
!  For each reflex vertex, resolve it with one or two separators
!  and update VCL, HVL, PVL, IANG.
!
  piptol = pi + tol
  v = 1

  do

    if ( nvert < v ) then
      exit
    end if

    if ( piptol < iang(v) ) then

      call resvrt ( v, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, &
        maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call insed2 ( v ,w1, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
        pvl, iang, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      if ( 0 < w2 ) then
        call insed2 ( v, w2, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
          pvl, iang, ierror )
      end if

      if ( ierror /= 0 ) then
        return
      end if

    end if

    v = v + 1

  end do

  return
end
