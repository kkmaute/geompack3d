subroutine mfdec2 ( hflag, umdf, kappa, angspc, angtol, dmin, nmin, ntrid, &
  nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl, &
  pvl, iang, ivrt, xivrt, widsq, edgval, vrtval, area, psi, iwk, wk, ierr )

!*****************************************************************************80
!
!! MFDEC2 further divides convex polygons to limit mesh function variation.
!
!  Discussion:
!
!    This routine further subdivides convex polygons so that the variation
!    of a heuristic or user-supplied mesh distribution function in
!    each polygon is limited.
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
!    Input, real ( kind = 8 ) KAPPA, the mesh smoothness parameter in interval
!    [0.0,1.0].
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in radians
!    used to determine extra points as possible endpoints of separators.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter in
!    radians used in accepting separators.
!
!    Input, real ( kind = 8 ) DMIN, the parameter used to determine if
!    variation of mdf in polygon is 'sufficiently high'.
!
!    Input, integer ( kind = 4 ) NMIN, the parameter used to determine if 'sufficiently
!    large' number of triangles in polygon.
!
!    Input, integer ( kind = 4 ) NTRID, the desired number of triangles in mesh.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or
!    positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions or
!    positions used in HVL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of polygon vertices or
!    positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM,
!    AREA, PSI arrays.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 3*NVRT + INT(2*PI/ANGSPC) where NVRT is maximum number of
!    vertices in a convex polygon of the (input) decomposition.
!
!    Input, integer ( kind = 4 ) MAXWK - maximum size available for WK array; should be about
!    NPOLG + 3*(NVRT + INT(2*PI/ANGSPC)) + 2.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) REGNUM(1:NPOLG), the region numbers.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Input, integer ( kind = 4 ) IVRT(1:NVERT), XIVRT(1:NPOLG+1), real ( kind = 8 )
!    WIDSQ(1:NPOLG), EDGVAL(1:NVERT), VRTVAL(1:NVC), arrays output from
!    routine DSMDF2;  if .NOT. HFLAG then only first two arrays exist.
!
!    Input/output, real ( kind = 8 ) AREA(1:NPOLG), the area of convex
!    polygons in decomposition.
!
!    Output, real ( kind = 8 ) PSI(1:NPOLG), the mean mdf values in the
!    convex polygons.
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
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) angsp2
  real    ( kind = 8 ) angspc
  real    ( kind = 8 ) angtol
  real    ( kind = 8 ) area(maxhv)
  real    ( kind = 8 ) areapg
  real    ( kind = 8 ) arearg
  real    ( kind = 8 ) c1
  real    ( kind = 8 ) c2
  real    ( kind = 8 ) cosalp
  real    ( kind = 8 ) ctrx
  real    ( kind = 8 ) ctry
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) dmin
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  integer ( kind = 4 ), parameter :: edgv = 4
  real    ( kind = 8 ) edgval(nvert)
  logical              hflag
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) indpvl
  real    ( kind = 8 ) intreg
  integer ( kind = 4 ) ivrt(nvert)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) kappa
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) maxn
  real    ( kind = 8 ) mdfint
  integer ( kind = 4 ) mdftr
  real    ( kind = 8 ) mean
  integer ( kind = 4 ) nev
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) np
  integer ( kind = 4 ) ntrid
  real    ( kind = 8 ) numer
  integer ( kind = 4 ) nvrt
  real    ( kind = 8 ) nwarea
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pi2
  integer ( kind = 4 ), parameter :: polg = 2
  real    ( kind = 8 ) psi(maxhv)
  integer ( kind = 4 ) pvl(4,maxpv)
  real    ( kind = 8 ) r
  integer ( kind = 4 ) regnum(maxhv)
  real    ( kind = 8 ) sinalp
  real    ( kind = 8 ) stdv
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sumx
  real    ( kind = 8 ) sumy
  real    ( kind = 8 ) theta1
  real    ( kind = 8 ) theta2
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), external :: umdf
  integer ( kind = 4 ) v
  real    ( kind = 8 ) vcl(2,maxvc)
  real    ( kind = 8 ) vrtval(nvc)
  integer ( kind = 4 ) w
  real    ( kind = 8 ) widsq(npolg)
  real    ( kind = 8 ) wk(maxwk)
  real    ( kind = 8 ) wsq
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xivrt(npolg+1)
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  integer ( kind = 4 ) yc
!
!  WK(1:NPOLG) is used for mdf standard deviation in polygons.
!  Compute AREARG = area of region and INTREG = estimated integral
!  of MDF2(X,Y) or UMDF(X,Y).
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  nvrt = 0

  do i = 1, npolg
    nvrt = max ( nvrt, xivrt(i+1) - xivrt(i) )
  end do

  if ( hflag .and. maxiw < 2 * nvrt ) then
    ierr = 6
    return
  else if ( maxwk < npolg + 3*nvrt + 2 ) then
    ierr = 7
    return
  end if

  listev = 1
  xc = npolg + 1
  yc = xc + nvrt + 1
  mdftr = yc + nvrt + 1
  arearg = 0.0D+00
  intreg = 0.0D+00
  nev = -1

  do i = 1, npolg

    if ( hflag ) then
      wsq = widsq(i)
      call prmdf2 ( i, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, iwk(listev) )
    end if

    if ( nev == 0 ) then

      psi(i) = 1.0D+00 / wsq
      wk(i) = 0.0D+00
      mdfint = psi(i) * area(i)

    else

      nvrt = xivrt(i+1) - xivrt(i)
      k = xivrt(i)
      sumx = 0.0D+00
      sumy = 0.0D+00

      do j = 0,nvrt-1
        l = ivrt(k)
        wk(xc+j) = vcl(1,l)
        wk(yc+j) = vcl(2,l)
        sumx = sumx + wk(xc+j)
        sumy = sumy + wk(yc+j)
        k = k + 1
      end do

      ctrx = sumx / real ( nvrt, kind = 8 )
      ctry = sumy / real ( nvrt, kind = 8 )

      do j = 0,nvrt-1
        wk(xc+j) = wk(xc+j) - ctrx
        wk(yc+j) = wk(yc+j) - ctry
      end do

      wk(xc+nvrt) = wk(xc)
      wk(yc+nvrt) = wk(yc)
      call intpg ( nvrt, wk(xc), wk(yc), ctrx, ctry, area(i), hflag, umdf, &
        wsq, nev, ifv, iwk(listev), ivrt, edgval, vrtval, vcl, mdfint, &
        psi(i), wk(i), wk(mdftr) )

    end if

    arearg = arearg + area(i)
    intreg = intreg + mdfint

  end do
!
!  If HFLAG, compute mean mdf values from KAPPA, etc.  Scale PSI(I)'s
!  so that integral in region is 1. Determine which polygons need to
!  be further subdivided (indicated by negative PSI(I) value).
!
  if ( hflag ) then
    c1 = ( 1.0D+00 - kappa ) / intreg
    c2 = kappa / arearg
  else
    c1 = 1.0D+00 / intreg
    c2 = 0.0D+00
  end if

  do i = 1, npolg
    psi(i) = psi(i) * c1 + c2
    if ( psi(i) * dmin < c1 * wk(i) ) then
      if ( nmin < ntrid * psi(i) * area(i) ) then
        psi(i) = -psi(i)
      end if
    end if
  end do
!
!  Further subdivide polygons for which DMIN < STDV/MEAN and
!  (estimated number of triangles) greater than NMIN.
!
  angsp2 = 2.0D+00 * angspc
  pi2 = 2.0D+00 * pi
  inc = int ( pi2 / angspc )
  nev = 0
  np = npolg
  xc = 1

  do i = 1, np

    if ( psi(i) < 0.0D+00 ) then

      if ( hflag ) then
        wsq = widsq(i)
        call prmdf2 ( i, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, &
          iwk(listev) )
      end if

      l = npolg + 1
      k = i

60    continue

      if ( npolg < k ) go to 130

70    continue

      if ( 0.0D+00 <= psi(k) ) go to 120
      nvrt = 0
      sumx = 0.0D+00
      sumy = 0.0D+00
      j = hvl(k)

      do

        nvrt = nvrt + 1
        m = pvl(loc,j)
        sumx = sumx + vcl(1,m)
        sumy = sumy + vcl(2,m)
        j = pvl(succ,j)

        if ( j == hvl(k) ) then
          exit
        end if

      end do

      ctrx = sumx/ real ( nvrt, kind = 8 )
      ctry = sumy/ real ( nvrt, kind = 8 )
      maxn = nvrt + inc

      if ( maxiw < nev + maxn + 1 ) then
        ierr = 6
        return
      else if ( maxwk < 3*maxn + 2 ) then
        ierr = 7
        return
      end if

      yc = xc + maxn + 1
      mdftr = yc + maxn + 1
      indpvl = listev + nev
      nvrt = 0
      m = pvl(loc,j)
      x1 = vcl(1,m) - ctrx
      y1 = vcl(2,m) - ctry
      wk(xc) = x1
      wk(yc) = y1
      theta1 = atan2(y1,x1)
      p = j
      iwk(indpvl) = j

90    continue

      j = pvl(succ,j)
      m = pvl(loc,j)
      x2 = vcl(1,m) - ctrx
      y2 = vcl(2,m) - ctry
      theta2 = atan2(y2,x2)
      if ( theta2 < theta1) theta2 = theta2 + pi2
      delta = theta2 - theta1

      if ( angsp2 <= delta ) then

        m = int(delta/angspc)
        delta = delta/ real ( m, kind = 8 )
        dx = x2 - x1
        dy = y2 - y1
        numer = x1*dy - y1*dx
        alpha = theta1

        do ii = 1,m-1
          alpha = alpha + delta
          cosalp = cos(alpha)
          sinalp = sin(alpha)
          r = numer / ( dy * cosalp - dx * sinalp )
          nvrt = nvrt + 1
          wk(xc+nvrt) = r*cosalp
          wk(yc+nvrt) = r*sinalp
          iwk(indpvl+nvrt) = -p
        end do

      end if

      nvrt = nvrt + 1
      wk(xc+nvrt) = x2
      wk(yc+nvrt) = y2
      x1 = x2
      y1 = y2
      theta1 = theta2
      p = j
      iwk(indpvl+nvrt) = j

      if ( j /= hvl(k)) then
        go to 90
      end if

      call intpg ( nvrt, wk(xc), wk(yc), ctrx, ctry, area(k), hflag, &
        umdf, wsq, nev, ifv, iwk(listev), ivrt, edgval, vrtval, &
        vcl, mdfint, mean, stdv, wk(mdftr) )

      psi(k) = mean*c1 + c2

      if ( psi(k) * dmin < c1*stdv ) then

        if ( nmin < ntrid*psi(k) * area(k) ) then

          call sepmdf(angtol,nvrt,wk(xc),wk(yc),area(k), &
            mean,wk(mdftr),iwk(indpvl),iang,i1,i2)

          if ( i1 < 0 ) then

            if ( maxwk < yc + 3*nvrt ) then
              ierr = 7
              return
            end if

            call sepshp(angtol,nvrt,wk(xc),wk(yc), &
              iwk(indpvl),iang,i1,i2,wk(yc+nvrt+1), ierr )

            if ( ierr /= 0 ) then
              return
            end if

          end if

          if ( i1 < 0 ) then
            ierr = 222
            return
          end if

          v = iwk(indpvl+i1)

          if ( v < 0 ) then
            call insvr2(wk(xc+i1)+ctrx,wk(yc+i1)+ctry,-v, &
              nvc,nvert,maxvc,maxpv,vcl,pvl,iang,v,ierr)
            if ( ierr /= 0 ) then
              return
            end if
          end if

          w = iwk(indpvl+i2)

          if ( w < 0 ) then
            call insvr2(wk(xc+i2)+ctrx,wk(yc+i2)+ctry,-w, &
              nvc,nvert,maxvc,maxpv,vcl,pvl,iang,w,ierr)
            if ( ierr /= 0 ) then
              return
            end if
          end if

          call insed2(v,w,npolg,nvert,maxhv,maxpv,vcl, &
            regnum,hvl,pvl,iang,ierr)

          if ( ierr /= 0 ) then
            return
          end if

          nvrt = 0
          j = hvl(k)

          do

            m = pvl(loc,j)
            wk(xc+nvrt) = vcl(1,m)
            wk(yc+nvrt) = vcl(2,m)
            nvrt = nvrt + 1
            j = pvl(succ,j)
            if ( j == hvl(k) ) then
              exit
            end if

          end do

          nwarea = areapg ( nvrt, wk(xc), wk(yc) ) * 0.5D+00
          area(npolg) = area(k) - nwarea
          area(k) = nwarea
          psi(k) = -psi(k)
          psi(npolg) = psi(k)

        end if

      end if

      go to 70

120   continue

      if ( k == i ) then
        k = l
      else
        k = k + 1
      end if

      go to 60

130   continue

    end if

  end do

  return
end
