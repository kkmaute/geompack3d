subroutine insvr2 ( xi, yi, wp, nvc, nvert, maxvc, maxpv, vcl, pvl, iang, &
  w, ierr )

!*****************************************************************************80
!
!! INSVR2 inserts a point into the vertex coordinate and polygon vertex lists.
!
!  Discussion:
!
!    This routine inserts the point (XI,YI) into the vertex coordinate list and
!    polygon vertex list data structures.
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
!    Input, real ( kind = 8 ) XI, YI, the coordinates of point to be inserted.
!
!    Input, integer ( kind = 4 ) WP, the index of vertex in PVL which is to be the
!    predecessor vertex of the inserted vertex.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), polygon vertex list.
!
!    Input/output, real ( kind = 8 ) IANG(1:NVERT), the polygon interior angles.
!
!    Output, integer ( kind = 4 ) W, the index of inserted vertex in PVL.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc

  integer ( kind = 4 ), parameter :: edgv = 4
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) w
  integer ( kind = 4 ) wp
  integer ( kind = 4 ) ws
  integer ( kind = 4 ) ww
  integer ( kind = 4 ) wwp
  integer ( kind = 4 ) wws
  real    ( kind = 8 ) xi
  real    ( kind = 8 ) yi

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( maxvc <= nvc ) then
    ierr = 3
    return
  else if ( maxpv < nvert+2 ) then
    ierr = 5
    return
  end if
!
!  Update linked list of vertices of polygon containing vertex WP.
!
  nvc = nvc + 1
  vcl(1,nvc) = xi
  vcl(2,nvc) = yi
  ws = pvl(succ,wp)
  nvert = nvert + 1
  w = nvert
  pvl(loc,w) = nvc
  pvl(polg,w) = pvl(polg,wp)
  pvl(succ,wp) = w
  pvl(succ,w) = ws
  iang(w) = pi
  pvl(edgv,w) = pvl(edgv,wp)
!
!  If edge containing (XI,YI) is shared by another polygon,
!  then also update linked list of vertices of that polygon.
!
  if ( 0 < pvl(edgv,wp) ) then
    wws = pvl(edgv,wp)
    wwp = pvl(succ,wws)
    nvert = nvert + 1
    ww = nvert
    pvl(loc,ww) = nvc
    pvl(polg,ww) = pvl(polg,wws)
    pvl(succ,wws) = ww
    pvl(succ,ww) = wwp
    iang(ww) = pi
    pvl(edgv,wp) = ww
    pvl(edgv,ww) = wp
    pvl(edgv,wws) = w
  end if

  return
end
