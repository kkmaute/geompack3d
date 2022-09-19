subroutine jnhole ( itophv, angspc, angtol, nvc, nvert, maxvc, maxpv, &
  maxiw, maxwk, vcl, hvl, pvl, iang, iwk, wk, ierr )

!*****************************************************************************80
!
!! JNHOLE joins a hole boundary to the boundary of the surrounding polygon.
!
!  Discussion:
!
!    This routine joins a hole boundary to the boundary of the polygon
!    containing the hole by finding a cut edge originating from the top
!    vertex of hole to a point on outer polygon boundary above it.
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
!    Input, integer ( kind = 4 ) ITOPHV, the index in PVL of top vertex of hole.
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter used
!    in controlling the vertices to be considered as an endpoint of a separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter used
!    in accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 3 times number of vertices in outer polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 5 times number of vertices in outer polygon.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:*), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) angle
  real    ( kind = 8 ) angspc
  real    ( kind = 8 ) angtol
  real    ( kind = 8 ) dy
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) hv
  integer ( kind = 4 ) hvl(*)
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ilft
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) irgt
  integer ( kind = 4 ) itophv
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivs
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lw
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  real    ( kind = 8 ) s
  real    ( kind = 8 ) slft
  real    ( kind = 8 ) srgt
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) succil
  integer ( kind = 4 ) succir
  real    ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vp
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) vs
  integer ( kind = 4 ) vv
  integer ( kind = 4 ) w
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) ww
  real    ( kind = 8 ) xint
  real    ( kind = 8 ) xlft
  real    ( kind = 8 ) xrgt
  real    ( kind = 8 ) xt
  real    ( kind = 8 ) xv
  real    ( kind = 8 ) xvs
  real    ( kind = 8 ) ylft
  real    ( kind = 8 ) yrgt
  real    ( kind = 8 ) yt
  real    ( kind = 8 ) ytmtol
  real    ( kind = 8 ) ytptol
  real    ( kind = 8 ) yv
  real    ( kind = 8 ) yvs

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( maxvc < nvc+3 ) then
    ierr = 3
    return
  else if ( maxpv < nvert+5 ) then
    ierr = 5
    return
  end if
!
!  Determine 'closest' vertices on outer boundary which are to the
!  left and right of top vertex of hole and on the horizontal line
!  through top vertex. The two closest vertices must be on edges
!  which intersect the horizontal line and are partially above the
!  line. Ties are broken (in the case of a vertex on a cut edge)
!  by choosing the vertex on the edge of maximum or minimum dx/dy
!  slope depending on whether the vertex is to the left or right
!  of top vertex, respectively.
!
  ipoly = pvl(polg,itophv)
  lv = pvl(loc,itophv)
  xt = vcl(1,lv)
  yt = vcl(2,lv)
  dy = 0.0D+00
  hv = hvl(ipoly)
  iv = hv
  yv = vcl(2,pvl(loc,iv))

  do

    iv = pvl(succ,iv)
    yvs = vcl(2,pvl(loc,iv))
    dy = max ( dy, abs ( yvs - yv ) )
    yv = yvs

    if ( iv == hv ) then
      exit
    end if

  end do

  ytmtol = yt - tol*dy
  ytptol = yt + tol*dy
  ilft = 0
  irgt = 0
  xlft = 0.0D+00
  xrgt = 0.0D+00
  hv = hvl(ipoly)
  iv = hv
  lv = pvl(loc,iv)
  xv = vcl(1,lv)
  yv = vcl(2,lv)

20 continue

  ivs = pvl(succ,iv)
  lv = pvl(loc,ivs)
  xvs = vcl(1,lv)
  yvs = vcl(2,lv)

  if ( yv <= ytptol .and. ytptol < yvs ) then

    if ( ytmtol <= yv ) then

      if ( xt < xv ) then

        if ( xv < xrgt .or. irgt == 0 ) then

          irgt = iv
          xrgt = xv
          yrgt = yv
          srgt = (xvs - xv) / (yvs - yv)

        else if ( xv == xrgt ) then

          s = (xvs - xv) / (yvs - yv)

          if ( s < srgt ) then
            irgt = iv
            yrgt = yv
            srgt = s
          end if

        end if

      end if

    else

      xint = (yt - yv) * (xvs - xv) / (yvs - yv) + xv

      if ( xt < xint ) then

        if ( xint < xrgt .or. irgt == 0 ) then
          irgt = iv
          xrgt = xint
          yrgt = yt
        end if

      end if

    end if

  else if ( ytptol < yv .and. yvs <= ytptol ) then

    if ( ytmtol <= yvs ) then

      if ( xvs < xt ) then

        if ( xlft < xvs .or. ilft == 0 ) then

          ilft = iv
          xlft = xvs
          ylft = yvs
          slft = ( xvs - xv ) / ( yvs - yv )

        else if ( xvs == xlft ) then

          s = ( xvs - xv ) / ( yvs - yv )

          if ( slft < s ) then
            ilft = iv
            ylft = yvs
            slft = s
          end if

        end if

      end if

    else

      xint = ( yt - yv ) * ( xvs - xv ) / ( yvs - yv ) + xv

      if ( xint < xt ) then

        if ( xlft < xint .or. ilft == 0 ) then
          ilft = iv
          xlft = xint
          ylft = yt
        end if

      end if

    end if

  end if

  iv = ivs
  xv = xvs
  yv = yvs

  if ( iv /= hv ) then
    go to 20
  end if

  if ( ilft == 0 .or. irgt == 0 ) then
    ierr = 218
    return
  end if
!
!  Temporarily modify PVL to pass the subregion 'above' top vertex
!  of hole to routine RESVRT. The top vertex is the reflex vertex
!  passed to RESVRT (in the temporary subregion, it has interior
!  angle PI). This causes one separator to be chosen by RESVRT
!  and its other endpoint is above the top vertex.
!
  succil = pvl(succ,ilft)
  succir = pvl(succ,irgt)
  vcl(1,nvc+2) = xlft
  vcl(2,nvc+2) = ylft
  vcl(1,nvc+3) = xrgt
  vcl(2,nvc+3) = yrgt
  vp = nvert + 3
  vr = nvert + 4
  vs = nvert + 5
  iang(vr) = angle ( xlft, ylft, xt, yt, xrgt, yrgt )

  if ( iang(vr) < pi - tol .or. pi + tol < iang(vr) ) then
    ierr = 219
    return
  end if

  pvl(loc,vp) = nvc + 2
  pvl(polg,vp) = ipoly
  pvl(succ,vp) = vr
  pvl(edgv,vp) = 0
  pvl(loc,vr) = pvl(loc,itophv)
  pvl(polg,vr) = ipoly
  pvl(succ,vr) = vs
  pvl(edgv,vr) = 0
  pvl(loc,vs) = nvc + 3
  pvl(polg,vs) = ipoly
  pvl(succ,vs) = succir
  pvl(edgv,vs) = pvl(edgv,irgt)
  pvl(succ,ilft) = vp
  lv = pvl(loc,ilft)
  iang(vp) = angle(vcl(1,lv),vcl(2,lv),xlft,ylft,xt,yt)
  lv = pvl(loc,succir)
  iang(vs) = angle(xt,yt,xrgt,yrgt,vcl(1,lv),vcl(2,lv))
  w = 0

  call resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw, maxwk, &
    vcl, pvl, iang, w, ww, iwk, wk, ierr )
!
!  Remove temporary modification to PVL. There are three cases
!  depending on where the endpoint of separator is located:
!  successor of closest vertex to the right of top vertex,
!  predecessor of closest vertex to the left of top vertex,
!  or neither of these.
!
  if ( pvl(succ,vs) == w ) then
    pvl(succ,ilft) = succil
    pvl(succ,irgt) = w
    pvl(edgv,irgt) = pvl(edgv,vs)
    if ( 0 < pvl(edgv,irgt) ) then
      pvl(edgv,pvl(edgv,irgt)) = irgt
    end if
  else if ( pvl(succ,ilft) == w ) then
    pvl(succ,w) = succil
  else
    pvl(succ,ilft) = succil
  end if

  if ( ierr /= 0 ) then
    return
  end if
!
!  Update PVL with cut edge, i.e. join linked lists of vertices
!  of the hole polygon and the outer boundary polygon into one
!  linked list of vertices by adding the cut edge from the top
!  vertex of hole to the vertex on the outer boundary.
!
  nvert = nvert + 2
  vv = nvert - 1
  ww = nvert
  lv = pvl(loc,itophv)
  lw = pvl(loc,w)
  pvl(loc,vv) = lv
  pvl(loc,ww) = lw
  pvl(polg,vv) = ipoly
  pvl(polg,ww) = ipoly
  pvl(succ,vv) = pvl(succ,itophv)
  pvl(succ,ww) = pvl(succ,w)
  pvl(succ,itophv) = ww
  pvl(succ,w) = vv
  pvl(edgv,vv) = pvl(edgv,itophv)
  pvl(edgv,ww) = pvl(edgv,w)
  pvl(edgv,itophv) = w
  pvl(edgv,w) = itophv

  if ( 0 < pvl(edgv,vv) ) then
    pvl(edgv,pvl(edgv,vv)) = vv
  end if

  if ( 0 < pvl(edgv,ww) ) then
    pvl(edgv,pvl(edgv,ww)) = ww
  end if

  l = pvl(loc,pvl(succ,vv))
  iang(vv) = angle(vcl(1,lw),vcl(2,lw),vcl(1,lv),vcl(2,lv), &
    vcl(1,l),vcl(2,l))
  iang(itophv) = iang(itophv) - iang(vv)
  l = pvl(loc,pvl(succ,ww))
  iang(ww) = angle(vcl(1,lv),vcl(2,lv),vcl(1,lw),vcl(2,lw), &
    vcl(1,l),vcl(2,l))
  iang(w) = iang(w) - iang(ww)

  if ( msglvl == 2 ) then
    write ( *,600) itophv,w,vcl(1,lv),vcl(2,lv), vcl(1,lw),vcl(2,lw)
  end if

  600 format (1x,2i7,4f15.7)

  return
end
