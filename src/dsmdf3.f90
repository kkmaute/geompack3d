subroutine dsmdf3 ( nvc, nface, nvert, npolh, maxiw, maxwk, vcl, facep, &
  nrml, fvl, eang, hfl, pfl, ivrt, xivrt, ifac, xifac, wid, facval, edgval, &
  vrtval, ncface, ncedge, iwk, wk, ierr )

!*****************************************************************************80
!
!! DSMDF3 sets up a mesh distribution function data structure in 3D.
!
!  Discussion:
!
!    This routine sets up a data structure for heuristic mesh distribution
!    function from convex polyhedron decomposition data structure.
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
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates in VCL array.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces or positions used in FACEP array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of positions used in FVL array.
!
!    Input, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions used in
!    HFL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    greater than or equal to MAX ( NFACE, NCFACE + 6*NCVERT ) where
!    NCFACE = maximum number of faces in a polyhedron,
!    NCVERT = 2 * maximum number edges in a polyhedron.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    greater than or equal to 3*NCFACE + NCVERT.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal vectors
!    for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Input, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices in
!    PFL for each polyhedron.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:*), the list of signed face indices for
!    each polyhedron.
!
!    Output, integer ( kind = 4 ) IVRT(1:NVERT), the indices of face vertices in VCL,
!    ordered by face.
!
!    Output, integer ( kind = 4 ) XIVRT(1:NFACE+1), the pointer to first vertex of
!    each face in IVRT; vertices of face K are IVRT(I) for I from XIVRT(K)
!    to XIVRT(K+1)-1.
!
!    Output, integer ( kind = 4 ) IFAC(1:*), the indices of polyhedron faces in
!    FACEP, ordered by polyhedron; same size as PFL.
!
!    Output, integer ( kind = 4 ) XIFAC(1:NPOLH+1), the pointer to first face of each
!    polyhedron in IFAC; faces of polyhedron K are IFAC(I) for I from
!    XIFAC(K) to XIFAC(K+1)-1.
!
!    Output, real ( kind = 8 ) WID(1:NPOLH), the width of convex polyhedra.
!
!    Output, real ( kind = 8 ) FACVAL(1:NFACE), the value associated with
!    each face of decomposition.
!
!    Output, real ( kind = 8 ) EDGVAL(1:NVERT), the value associated with
!    each edge of decomposition.
!
!    Output, real ( kind = 8 ) VRTVAL(1:NVC), the value associated with
!    each vertex of decomposition.
!
!    Output, integer ( kind = 4 ) NCFACE, the maximum number of faces in a polyhedron.
!
!    Output, integer ( kind = 4 ) NCEDGE, the maximum number of edges in a polyhedron.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npolh
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert

  integer ( kind = 4 ) ccw
  integer ( kind = 4 ) ceang
  integer ( kind = 4 ) cfvl
  integer ( kind = 4 ) chvl
  integer ( kind = 4 ) cnrml
  real    ( kind = 8 ) cxy
  real    ( kind = 8 ) cyz
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) dirsq
  real    ( kind = 8 ) dir1(3)
  real    ( kind = 8 ) dir1sq
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) eang(nvert)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) edgval(nvert)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,nface)
  integer ( kind = 4 ), parameter :: facn = 2
  real    ( kind = 8 ) facval(nface)
  real    ( kind = 8 ) fmax
  integer ( kind = 4 ) fvl(6,nvert)
  integer ( kind = 4 ) hfl(npolh)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac(*)
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) ivrt(nvert)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) li
  integer ( kind = 4 ) lj
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) ncedge
  integer ( kind = 4 ) ncface
  integer ( kind = 4 ) ncvert
  real    ( kind = 8 ) nrml(3,nface)
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,*)
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) r21
  real    ( kind = 8 ) r22
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sxy
  real    ( kind = 8 ) syz
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,nvc)
  real    ( kind = 8 ) vrtval(nvc)
  real    ( kind = 8 ) wid(npolh)
  real    ( kind = 8 ) widsq
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xifac(npolh+1)
  integer ( kind = 4 ) xivrt(nface+1)
  integer ( kind = 4 ) yc
!
!  Compute width of polyhedra.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( maxiw < nface ) then
    ierr = 6
    return
  end if

  ncface = 0
  ncvert = 0
  nvrt = 0

  do f = 1, nface

    k = 0
    i = facep(1,f)

    do

      k = k + 1
      i = fvl(succ,i)

      if ( i == facep(1,f) ) then
        exit
      end if

    end do

    iwk(f) = k
    nvrt = max ( nvrt, k )

  end do

  do p = 1, npolh

    j = 0
    k = 0
    i = hfl(p)

    do

      f = abs ( pfl(1,i) )
      j = j + 1
      k = k + iwk(f)
      i = pfl(2,i)

      if ( i == hfl(p) ) then
        exit
      end if

    end do

    ncface = max ( ncface, j )
    ncvert = max ( ncvert, k )

  end do

  ncedge = ncvert / 2
  chvl = 1
  cfvl = chvl + ncface
  irem = cfvl + 5 * ncvert
  cnrml = 1
  ceang = cnrml + 3 * ncface

  if ( maxiw < irem + ncvert - 1 ) then
    ierr = 6
    return
  else if ( maxwk < ceang + ncvert - 1 ) then
    ierr = 7
    return
  end if

  do i = 1, npolh

    call dsconv(i,hfl(i),facep,nrml,fvl,eang,pfl,ncface,ncvert, &
      iwk(chvl),wk(cnrml),iwk(cfvl),wk(ceang))

    irem = cfvl + 5*ncvert

    call rmcpfc(ncface,ncvert,iwk(chvl),wk(cnrml),iwk(cfvl), &
      wk(ceang),iwk(irem))

    call rmcled(ncface,ncvert,iwk(chvl),iwk(cfvl))

    irem = cfvl + 5*ncvert

    call width3(ncface,vcl,iwk(chvl),wk(cnrml),iwk(cfvl), &
      maxiw-irem+1,j,k,wid(i),iwk(irem), ierr )

    if ( ierr /= 0 ) then
      return
    end if

  end do
!
!  Set up FACVAL array.
!  For each face, rotate normal vector to (0,0,1). Rotation matrix is
!    [ CXY     -SXY     0   ]
!    [ CYZ*SXY CYZ*CXY -SYZ ]
!    [ SYZ*SXY SYZ*CXY  CYZ ]
!  Rotate face vertices with int angles < PI before call to WIDTH2.
!
  fmax = 0.0D+00
  xc = 1
  yc = xc + nvrt

  do f = 1, nface

    if ( abs(nrml(1,f)) <= tol ) then
      leng = nrml(2,f)
      cxy = 1.0D+00
      sxy = 0.0D+00
    else
      leng = sqrt(nrml(1,f)**2 + nrml(2,f)**2)
      cxy = nrml(2,f)/leng
      sxy = nrml(1,f)/leng
    end if

    cyz = nrml(3,f)
    syz = leng
    r21 = cyz*sxy
    r22 = cyz*cxy

    if ( 0 < facep(2,f) ) then
      ccw = succ
    else
      ccw = pred
    end if

    nvrt = 0
    i = facep(1,f)
    li = fvl(loc,i)
    lj = fvl(loc,fvl(7-ccw,i))
    dir(1) = vcl(1,li) - vcl(1,lj)
    dir(2) = vcl(2,li) - vcl(2,lj)
    dir(3) = vcl(3,li) - vcl(3,lj)
    dirsq = dir(1)**2 + dir(2)**2 + dir(3)**2

60  continue

    j = fvl(ccw,i)
    lj = fvl(loc,j)
    dir1(1) = vcl(1,lj) - vcl(1,li)
    dir1(2) = vcl(2,lj) - vcl(2,li)
    dir1(3) = vcl(3,lj) - vcl(3,li)
    dir1sq = dir1(1)**2 + dir1(2)**2 + dir1(3)**2
    dotp = -(dir(1)*dir1(1) + dir(2)*dir1(2) + dir(3)*dir1(3))/ &
      sqrt(dirsq*dir1sq)

    if ( -1.0D+00 + tol < dotp ) then
      wk(xc+nvrt) = cxy*vcl(1,li) - sxy*vcl(2,li)
      wk(yc+nvrt) = r21*vcl(1,li)+r22*vcl(2,li)-syz*vcl(3,li)
      nvrt = nvrt + 1
    end if

    i = j
    li = lj
    dir(1) = dir1(1)
    dir(2) = dir1(2)
    dir(3) = dir1(3)
    dirsq = dir1sq

    if ( i /= facep(1,f) ) then
      go to 60
    end if

    call width2(nvrt,wk(xc),wk(yc),i,j,widsq, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    facval(f) = min ( sqrt ( widsq ), wid(abs(facep(2,f))) )

    if ( facep(3,f) /= 0 ) then
      facval(f) = min ( facval(f), wid(abs(facep(3,f))) )
    end if

    fmax = max ( fmax, facval(f) )

  end do
!
!  Set up IFAC, XIFAC, IVRT, XIVRT arrays.
!
  k = 1

  do p = 1, npolh

    xifac(p) = k
    i = hfl(p)

    do

      ifac(k) = pfl(1,i)
      i = pfl(2,i)
      k = k + 1

      if ( i == hfl(p) ) then
        exit
      end if

    end do

  end do

  xifac(npolh+1) = k

  k = 1

  do f = 1, nface

    xivrt(f) = k
    i = facep(1,f)

    do

      ivrt(k) = fvl(loc,i)
      fvl(pred,i) = k
      i = fvl(succ,i)
      k = k + 1

      if ( i == facep(1,f) ) then
        exit
      end if

    end do

  end do

  xivrt(nface+1) = k
!
!  Set up EDGVAL, VRTVAL arrays and reset FVL(PRED,*).
!
  edgval(1:nvert) = 0.0D+00

  vrtval(1:nvc) = fmax

  do i = 1, nvert

    if ( 0.0D+00 < edgval(fvl(pred,i)) ) then
      cycle
    end if

    li = fvl(loc,i)
    lj = fvl(loc,fvl(succ,i))
    leng = sqrt((vcl(1,lj) - vcl(1,li))**2 + (vcl(2,lj) - &
      vcl(2,li))**2 + (vcl(3,lj) - vcl(3,li))**2)
    k = i
    j = i

140 continue

    f = fvl(facn,j)
    leng = min ( leng, facval(f) )

    if ( fvl(edga,j) == 0 ) then

      k = j
      j = fvl(edgc,i)

      do while ( j /= 0 )
        f = fvl(facn,j)
        leng = min ( leng, facval(f) )
        j = fvl(edgc,j)
      end do

    else if ( fvl(edga,j) /= i ) then

      j = fvl(edga,j)
      go to 140

    end if

    j = k

160 continue

    edgval(fvl(pred,j)) = leng
    j = fvl(edgc,j)

    if ( j /= k .and. j /= 0 ) then
      go to 160
    end if

    vrtval(li) = min ( vrtval(li), leng )
    vrtval(lj) = min ( vrtval(lj), leng )

  end do

  do f = 1, nface

    i = facep(1,f)

    do

      j = fvl(succ,i)
      fvl(pred,j) = i
      i = j

      if ( i == facep(1,f) ) then
        exit
      end if

    end do

  end do

  return
end
