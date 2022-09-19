subroutine insed2 ( v, w, npolg, nvert, maxhv, maxpv, vcl, regnum, hvl, &
  pvl, iang, ierr )

!*****************************************************************************80
!
!! INSED2 inserts an edge into the head and polygon vertex lists.
!
!  Discussion:
!
!    This routine inserts an edge joining vertices V, W into head vertex
!    list and polygon vertex list data structures.
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
!    Input, integer ( kind = 4 ) V, W, indices in PVL of vertices which are the endpoints
!    of an edge to be added to decomposition.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of positions used in HVL array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), polygon vertex list.
!
!    Input/output, real ( kind = 8 ) IANG(1:NVERT), interior angles.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxpv

  real    ( kind = 8 ) angle
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lw
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) v
  real    ( kind = 8 ) vcl(2,*)
  integer ( kind = 4 ) vv
  integer ( kind = 4 ) w
  integer ( kind = 4 ) ww

  ierr = 0

  if ( maxhv <= npolg ) then
    ierr = 4
    return
  else if ( maxpv < nvert+2 ) then
    ierr = 5
    return
  end if
!
!  Split linked list of vertices of the polygon containing vertices
!  V and W into two linked list of vertices of polygons with common
!  edge joining V and W.
!
  nvert = nvert + 2
  vv = nvert - 1
  ww = nvert
  lv = pvl(loc,v)
  lw = pvl(loc,w)
  pvl(loc,vv) = lv
  pvl(loc,ww) = lw
  pvl(polg,ww) = pvl(polg,v)
  pvl(succ,vv) = pvl(succ,v)
  pvl(succ,ww) = pvl(succ,w)
  pvl(succ,v) = ww
  pvl(succ,w) = vv
  pvl(edgv,vv) = pvl(edgv,v)
  pvl(edgv,ww) = pvl(edgv,w)
  pvl(edgv,v) = w
  pvl(edgv,w) = v

  if ( 0 < pvl(edgv,vv) ) then
    pvl(edgv,pvl(edgv,vv)) = vv
  end if

  if ( 0 < pvl(edgv,ww) ) then
    pvl(edgv,pvl(edgv,ww)) = ww
  end if

  l = pvl(loc,pvl(succ,vv))
  iang(vv) = angle(vcl(1,lw),vcl(2,lw),vcl(1,lv),vcl(2,lv),vcl(1,l),vcl(2,l))
  iang(v) = iang(v) - iang(vv)
  l = pvl(loc,pvl(succ,ww))
  iang(ww) = angle(vcl(1,lv),vcl(2,lv),vcl(1,lw),vcl(2,lw),vcl(1,l),vcl(2,l))
  iang(w) = iang(w) - iang(ww)
  npolg = npolg + 1
  i = vv

  do

    pvl(polg,i) = npolg
    i = pvl(succ,i)

    if ( i == vv ) then
      exit
    end if

  end do

  hvl(pvl(polg,v)) = v
  hvl(npolg) = vv
  regnum(npolg) = regnum(pvl(polg,v))

  if ( msglvl == 2 ) then
    write ( *, '(2i7,4f15.7)' ) v,w,vcl(1,lv),vcl(2,lv), vcl(1,lw),vcl(2,lw)
  end if

  return
end
