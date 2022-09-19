subroutine cvdecf ( aspc2d, atol2d, nvc, nface, nvert, npf, maxvc, maxfp, &
  maxfv, maxpf, maxiw, maxwk, vcl, facep, factyp, nrml, fvl, eang, hfl, pfl, &
  iwk, wk, ierr )

!*****************************************************************************80
!
!! CVDECF updates a polyhedral decomposition.
!
!  Discussion:
!
!    This routine updates a polyhedral decomposition data structure for
!    a polyhedral region by decomposing each face of the decomposition into
!    convex subpolygons using separators from reflex vertices.
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
!    Input, real ( kind = 8 ) ASPC2D, the angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of a
!    separator.
!
!    Input, real ( kind = 8 ) ATOL2D, the angle tolerance parameter in
!    radians used in accepting separator to resolve a reflex vertex on a face.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or
!    positions used in VCL.
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces or positions used
!    in FACEP array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input/output, integer ( kind = 4 ) NPF, the number of positions used in PFL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXFP, the maximum size available for FACEP, FACTYP,
!    NRML arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXPF, the maximum size available for PFL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 6*NV + 24*NRV where NV is max number of vertices and NRV
!    is maximum number of reflex vertices in a nonconvex face.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 8*NV + 24*NRV.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list:
!    row 1 is head pointer, rows 2 and 3 are signed polyhedron indices.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types: useful for
!    specifying types of boundary faces; entries must be greater than or
!    equal to 0; any new interior faces (not part of previous face) has
!    face type set to 0.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal
!    vectors for faces; outward normal corresponds to counterclockwise
!    traversal of face from polyhedron with index |FACEP(2,F)|.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list; see
!    routine DSPHDC.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the angles at edges
!    common to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J),
!    determined by EDGC field.
!
!    Input, integer ( kind = 4 ) HFL(1:*), the head pointer to face indices in PFL for
!    each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices
!    for each polyhedron; row 2 used for link.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  integer ( kind = 4 ) a
  real    ( kind = 8 ) angle
  real    ( kind = 8 ) aspc2d
  real    ( kind = 8 ) atol2d
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) cp
  real    ( kind = 8 ) cxy
  real    ( kind = 8 ) cyz
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edgv
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) hfl(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iang
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) js
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) link
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) locfv
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) ls
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npf
  real    ( kind = 8 ) nrml(3,maxfp)
  integer ( kind = 4 ) nrv
  real    ( kind = 8 ) ntol
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) pfl(2,maxpf)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) piptol
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) pt(2,2)
  real    ( kind = 8 ) r21
  real    ( kind = 8 ) r22
  real    ( kind = 8 ) r31
  real    ( kind = 8 ) r32
  integer ( kind = 4 ) s
  integer ( kind = 4 ) size
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sxy
  real    ( kind = 8 ) syz
  integer ( kind = 4 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) u(3)
  real    ( kind = 8 ) v(3)
  real    ( kind = 8 ) vcl(3,maxvc)
  integer ( kind = 4 ) wrem
  integer ( kind = 4 ) w(2)
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  real    ( kind = 8 ) zr

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  piptol = pi + tol
  nf = nface

  do i = 1, nf
!
!  Determine number of reflex vertices of face I.
!
    nrv = 0
    nv = 0
    k = 1
    if ( abs(nrml(1,i)) < abs(nrml(2,i))) k = 2
    if ( abs(nrml(k,i)) < abs(nrml(3,i))) k = 3

    if ( 0 < facep(2,i) ) then
      ccw = succ
    else
      ccw = pred
    end if

    j = facep(1,i)
    l = fvl(loc,j)
    lp = fvl(loc,fvl(7-ccw,j))
    u(1) = vcl(1,l) - vcl(1,lp)
    u(2) = vcl(2,l) - vcl(2,lp)
    u(3) = vcl(3,l) - vcl(3,lp)

10  continue

    nv = nv + 1
    js = fvl(ccw,j)
    ls = fvl(loc,js)
    v(1) = vcl(1,ls) - vcl(1,l)
    v(2) = vcl(2,ls) - vcl(2,l)
    v(3) = vcl(3,ls) - vcl(3,l)
    ntol = tol * max ( abs(u(1)), abs(u(2)), abs(u(3)), abs(v(1) ), &
      abs(v(2)), abs(v(3)) )

    if ( k == 1 ) then
      cp = u(2)*v(3) - u(3)*v(2)
    else if ( k == 2 ) then
      cp = u(3)*v(1) - u(1)*v(3)
    else
      cp = u(1)*v(2) - u(2)*v(1)
    end if

    if ( ntol < abs(cp) .and. cp*nrml(k,i) < 0.0D+00 ) then
      nrv = nrv + 1
    end if

    j = js

    if ( j /= facep(1,i) ) then
      l = ls
      u(1) = v(1)
      u(2) = v(2)
      u(3) = v(3)
      go to 10
    end if

    if ( nrv == 0 ) then
      go to 70
    end if
!
!  Set up 2D data structure. Rotate normal vector of face to
!  (0,0,1). Rotation matrix applied to face vertices is
!    [ CXY     -SXY     0   ]
!    [ CYZ*SXY CYZ*CXY -SYZ ]
!    [ SYZ*SXY SYZ*CXY  CYZ ]
!
    size = nv + 8*nrv
    locfv = 1
    link = locfv + size
    edgv = link + size
    irem = edgv + size
    x = 1
    y = x + size
    iang = y + size
    wrem = iang + size

    if ( maxiw < irem ) then
      ierr = 6
      return
    else if ( maxwk < wrem ) then
      ierr = 7
      return
    end if

    if ( abs(nrml(1,i)) <= tol ) then
      leng = nrml(2,i)
      cxy = 1.0D+00
      sxy = 0.0D+00
    else
      leng = sqrt(nrml(1,i)**2 + nrml(2,i)**2)
      cxy = nrml(2,i)/leng
      sxy = nrml(1,i)/leng
    end if

    cyz = nrml(3,i)
    syz = leng
    r21 = cyz*sxy
    r22 = cyz*cxy
    r31 = nrml(1,i)
    r32 = nrml(2,i)
    zr = r31*vcl(1,ls) + r32*vcl(2,ls) + cyz*vcl(3,ls)
    j = facep(1,i)

    do k = 0, nv-1
      l = fvl(loc,j)
      wk(x+k) = cxy*vcl(1,l) - sxy*vcl(2,l)
      wk(y+k) = r21*vcl(1,l) + r22*vcl(2,l) - syz*vcl(3,l)
      iwk(locfv+k) = j
      iwk(link+k) = k + 2
      iwk(edgv+k) = 0
      j = fvl(ccw,j)
    end do

    iwk(link+nv-1) = 1

    do k = 1, nv-2
      wk(iang+k) = angle(wk(x+k-1),wk(y+k-1), wk(x+k),wk(y+k), &
        wk(x+k+1),wk(y+k+1))
    end do

    wk(iang) = angle(wk(x+nv-1),wk(y+nv-1), wk(x),wk(y), &
      wk(x+1),wk(y+1))

    wk(iang+nv-1) = angle(wk(x+nv-2),wk(y+nv-2), wk(x+nv-1), &
      wk(y+nv-1), wk(x),wk(y))
!
!  Resolve reflex vertices.
!
    j = 1

40  continue

    if ( nv < j ) then
      go to 70
    end if

    if ( piptol < wk(iang+j-1) ) then

      call resvrf(j,aspc2d,atol2d,maxiw-irem+1,maxwk-wrem+1, &
        wk(x),wk(y),wk(iang),iwk(link),w(1),w(2),pt,pt(1,2), &
        iwk(irem),wk(wrem), ierr )

      if ( ierr /= 0 ) then
        return
      end if

      if ( w(2) == 0 ) then
        l = 1
      else
        l = 2
      end if

      do k = l, 1, -1

        s = -w(k)

        if ( 0 < s ) then

          wk(x+nv) = pt(1,k)
          wk(y+nv) = pt(2,k)
          wk(iang+nv) = pi
          iwk(link+nv) = iwk(link+s-1)
          nv = nv + 1
          iwk(link+s-1) = nv

          if ( maxvc <= nvc ) then
            ierr = 14
            return
          end if

          vcl(1,nvc+1) = cxy*pt(1,k) + r21*pt(2,k) + r31*zr
          vcl(2,nvc+1) = r22*pt(2,k) - sxy*pt(1,k) + r32*zr
          vcl(3,nvc+1) = cyz*zr - syz*pt(2,k)
          a = iwk(locfv+s-1)
          if ( ccw == pred) a = fvl(pred,a)
          call insvr3(a,nvc,nvert,maxfv,vcl,fvl,eang,ierr)

          if ( ierr /= 0 ) then
            return
          end if

          iwk(locfv+nv-1) = fvl(succ,a)
          t = iwk(edgv+s-1)

          if ( t == 0 ) then
            iwk(edgv+nv-1) = 0
            w(k) = nv
          else
            wk(x+nv) = wk(x+nv-1)
            wk(y+nv) = wk(y+nv-1)
            wk(iang+nv) = pi
            iwk(link+nv) = iwk(link+t-1)
            nv = nv + 1
            iwk(link+t-1) = nv

            if ( iwk(locfv+nv-2) == nvert ) then
              iwk(locfv+nv-1) = nvert - 1
            else
              iwk(locfv+nv-1) = nvert
            end if

            iwk(edgv+s-1) = nv
            iwk(edgv+nv-2) = t
            iwk(edgv+t-1) = nv - 1
            iwk(edgv+nv-1) = s
            w(k) = nv - 1
          end if

        end if

      end do

      do k = 1, l

        s = w(k)
        wk(x+nv) = wk(x+j-1)
        wk(y+nv) = wk(y+j-1)
        iwk(link+nv) = iwk(link+j-1)
        nv = nv + 1
        wk(x+nv) = wk(x+s-1)
        wk(y+nv) = wk(y+s-1)
        iwk(link+nv) = iwk(link+s-1)
        nv = nv + 1
        iwk(link+j-1) = nv
        iwk(link+s-1) = nv - 1
        a = iwk(edgv+j-1)
        iwk(edgv+nv-2) = a
        iwk(edgv+j-1) = s

        if ( 0 < a ) then
          iwk(edgv+a-1) = nv - 1
        end if

        a = iwk(edgv+s-1)
        iwk(edgv+nv-1) = a
        iwk(edgv+s-1) = j
        if ( 0 < a ) then
          iwk(edgv+a-1) = nv
        end if
        a = iwk(locfv+j-1)
        call insed3(a,iwk(locfv+s-1),nface,nvert,npf,maxfp, &
          maxfv,maxpf,facep,factyp,nrml,fvl,eang,hfl,pfl,ierr)

        if ( ierr /= 0 ) then
          return
        end if

        if ( ccw == succ ) then
          iwk(locfv+nv-2) = iwk(locfv+j-1)
          iwk(locfv+nv-1) = iwk(locfv+s-1)
          iwk(locfv+j-1) = nvert - 1
          iwk(locfv+s-1) = nvert
        else
          iwk(locfv+nv-2) = nvert - 1
          iwk(locfv+nv-1) = nvert
        end if

        t = iwk(link+nv-2)
        wk(iang+nv-2) = angle(wk(x+s-1),wk(y+s-1), &
                 wk(x+nv-2),wk(y+nv-2), wk(x+t-1),wk(y+t-1))
        wk(iang+j-1) = wk(iang+j-1) - wk(iang+nv-2)
        t = iwk(link+nv-1)
        wk(iang+nv-1) = angle(wk(x+j-1),wk(y+j-1), &
                 wk(x+nv-1),wk(y+nv-1), wk(x+t-1),wk(y+t-1))
        wk(iang+s-1) = wk(iang+s-1) - wk(iang+nv-1)

      end do

    end if

    j = j + 1
    go to 40

70  continue

  end do

  return
end
