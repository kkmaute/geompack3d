subroutine eqdis2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid, &
  nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl, &
  pvl, iang, area, psi, h, iwk, wk, ierr )

!*****************************************************************************80
!
!! EQDIS2 subdivides convex polygons for equidistribution.
!
!  Discussion:
!
!    This routine further subdivides convex polygons so that an approximate
!    equidistributing triangular mesh can be constructed with
!    respect to a heuristic or user-supplied mesh distribution
!    function, and determines triangle sizes for each polygon of
!    the decomposition.
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
!    Input, logical HFLAG, TRUE if heuristic mdf, FALSE if user-supplied mdf.
!
!    Input, external real ( kind = 8 ) UMDF(X,Y), the user-supplied mdf
!    with d.p arguments.
!
!    Input, real ( kind = 8 ) KAPPA, the mesh smoothness parameter in
!    interval [0.0,1.0].
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians
!    used to determine extra points as possible endpoints of separators.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in radians
!    used in accepting separators.
!
!    Input, real ( kind = 8 ) DMIN, the parameter used to determine if
!    variation of mdf in polygon is 'sufficiently high'.
!
!    Input, integer ( kind = 4 ) NMIN, the parameter used to determine if 'sufficiently
!    large' number of triangles in polygon.
!
!    Input, integer ( kind = 4 ) NTRID, the desired number of triangles in mesh.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or
!    positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array, should
!    be greater than or equal to the number of vertex coordinates required
!    for decomposition (approximately NVC + 2*NS where NS is expected number
!    of new separators).
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM,
!    AREA, PSI, H arrays; should be greater than or equal to the number
!    of polygons required for decomposition (approximately NPOLG + NS).
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays;
!    should be greater than or equal to the number of polygon vertices
!    required for decomposition (approximately NVERT + 5*NS).
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    greater than or equal to
!    MAX ( 2*NP, NVERT + NPOLG + 3*NVRT + INT(2*PI/ANGSPC) )
!    where NVRT is maximum number of vertices in a convex
!    polygon of the (input) decomposition, NP is expected
!    value of NPOLG on output.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    greater than or equal to
!    NVC + NVERT + 2*NPOLG + 3*(NVRT + INT(2*PI/ANGSPC)).
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate
!    list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), IANG(1:NVERT), the polygon
!    vertex list and interior angles; see routine DSPGDC for more details.
!
!    [Note: The data structures should be as output from routine CVDEC2.]
!
!    Output, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons
!    in decomposition.
!
!    Output, real ( kind = 8 ) PSI(1:NPOLG), the smoothed mean mdf values
!    in the convex polygons.
!
!    Output, real ( kind = 8 ) H(1:NPOLG), the triangle size for convex
!    polygons.
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
  real    ( kind = 8 ) area(maxhv)
  real    ( kind = 8 ) dmin
  integer ( kind = 4 ) edgval
  real    ( kind = 8 ) h(maxhv)
  logical              hflag
  integer ( kind = 4 ) hvl(maxhv)
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(maxiw)
  real    ( kind = 8 ) kappa
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) ntrid
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real    ( kind = 8 ) psi(maxhv)
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  real    ( kind = 8 ), external :: umdf
  real    ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vrtval
  integer ( kind = 4 ) widsq
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xivrt

  ierr = 0
  ivrt = 1
  xivrt = ivrt + nvert
  m = xivrt + npolg

  if ( maxiw < m ) then
    ierr = 6
    return
  end if

  widsq = 1

  if ( hflag ) then
    edgval = widsq + npolg
    vrtval = edgval + nvert
    n = npolg + nvert + nvc
    if ( maxwk < n ) then
      ierr = 7
      return
    end if
  else
    edgval = 1
    vrtval = 1
    n = 0
  end if

  call dsmdf2(hflag,nvc,npolg,maxwk-n,vcl,hvl,pvl,iang,iwk(ivrt), &
    iwk(xivrt),wk(widsq),wk(edgval),wk(vrtval),area,wk(n+1),ierr)

  if ( ierr /= 0 ) then
    return
  end if

  call mfdec2(hflag,umdf,kappa,angspc,angtol,dmin,nmin,ntrid,nvc, &
    npolg,nvert,maxvc,maxhv,maxpv,maxiw-m,maxwk-n,vcl,regnum,hvl, &
    pvl,iang,iwk(ivrt),iwk(xivrt),wk(widsq),wk(edgval),wk(vrtval), &
    area,psi,iwk(m+1),wk(n+1), ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( maxiw < 2 * npolg ) then
    ierr = 6
    return
  end if

  call trisiz(ntrid,npolg,hvl,pvl,area,psi,h,iwk,iwk(npolg+1))

  return
end
