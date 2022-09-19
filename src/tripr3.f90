subroutine tripr3 ( h, shrf, crit, nvc, nface, nvert, npolh, maxvc, maxbt, &
  fc_max, maxht, maxvm, maxiw, maxwk, vcl, facep, nrml, fvl, eang, hfl, pfl, &
  edst, edno, fcst, btst, btl, fc, ht, vm, xfhv, ntrif, ntetra, iwk, wk, ierr )

!*****************************************************************************80
!
!! TRIPR3 generates mesh vertices in a decomposed polygonal region.
!
!  Discussion:
!
!    This routine generates and triangulates mesh vertices in a polyhedral
!    region which has been decomposed into convex polyhedra with
!    convex faces.
!
!    In the P-th convex polyhedron, the mesh vertices are
!    generated on a quasi-uniform grid of spacing H(P). The initial
!    triangulation in each polyhedron is boundary-constrained based
!    on local empty circumsphere criterion (often Delaunay), and
!    this can be improved by local transformations based on max-min
!    solid angle or ratio of inradius to circumradius criterion.
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
!    Input, real ( kind = 8 ) H(1:NPOLH), the mesh spacings for convex
!    polyhedra in polyhedral decomposition data structure.
!
!    Input, real ( kind = 8 ) SHRF, the factor for shrinking polyhedron;
!    shrink distance is SHRF*H(P); a value of about 0.8 for SHRF is
!    recommended.
!
!    Input, integer ( kind = 4 ) CRIT, the code for improving b-c triangulations:
!    1 for (local max-min) solid angle criterion,
!    2 for radius ratio criterion,
!    3 for mean ratio crit, else no improvement.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces or positions used in FACEP array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of positions used in FVL, EANG array.
!
!    Input, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions used in
!    HFL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array; should
!    be large enough to hold all mesh vertices.
!
!    Input, integer ( kind = 4 ) MAXBT, the maximum size available for BTL array; should be
!    at least number of triangles generated on boundary faces.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array; should be
!    at least the sum of number of faces in triangulation of each polyhedron.
!
!    Input, integer ( kind = 4 ) MAXHT, the maximum size available for HT array; should be
!    at least 3/2 * MAXVM.
!
!    Input, integer ( kind = 4 ) MAXVM, the maximum size available for VM array; should be
!    at least the sum of number of mesh vertices in each polyhedron.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    at least NVCF + 2*NCFACE + 17*NCVERT where NVCF is number of mesh
!    vertices on all faces of decomposition (excluding those
!    in interior of polyhedra), NCFACE is max number of faces
!    in a polyhedron, NCVERT is 2 * max number of edges in a polyhedron.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    at least 3*NCFACE + 8*NCVERT.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list: row 1 is
!    head pointer, rows 2 and 3 are signed polyhedron indices.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal vectors
!    for faces; outward normal corresponds to counterclockwise traversal of
!    face from polyhedron with index |FACEP(2,F)|.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list; see routine DSPHDC.
!
!    Input, real ( kind = 8 ) EANG(1:NVERT), the angles at edges common
!    to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J), determined
!    by EDGC field.
!
!    Input, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices in PFL
!    for each polyhedron.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:*), the list of signed face indices for
!    each polyhedron; row 2 used for link.
!
!    Output, integer ( kind = 4 ) EDST(1:NVERT), the start location in VCL for mesh
!    vertices on each edge in FVL if there are any, else 0.
!
!    Output, integer ( kind = 4 ) EDNO(1:NVERT), the number of mesh vertices on interior
!    of each edge in FVL; entry is negated if mesh vertices are listed in
!    backward order (according to SUCC) in VCL.
!
!    Output, integer ( kind = 4 ) FCST(1:NFACE+1), the start location in VCL for mesh
!    vertices in interior of each face; last entry indicates where mesh
!    vertices in interior of polyhedra start.
!
!    Output, integer ( kind = 4 ) BTST(1:NFACE+1), the start location in BTL for triangles
!    on each face of FACEP; last entry is (total number of triangles + 1).
!
!    Output, integer ( kind = 4 ) BTL(1:3,1:NBT), the boundary triangle list for triangles
!    generated on all faces of decomposition; NBT = BTST(NFACE+1) - 1;
!    entries are indices of VCL.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), HT(1:MAXHT), VM(1:MAXVM), the face
!    record array, hash table, and vertex mapping array for triangulations
!    in each convex polyhedron (see routine BCDTRI).
!
!    Output, integer ( kind = 4 ) XFHV(1:3,1:NPOLH+1), the indicates usage of FC, HT,
!    VM arrays for polyhedron P; space used in FC for triangulation in P is
!    from FC(*,J:K) where J=XFHV(1,P) and K=XBFHV(1,P+1)-1; space
!    used for HT, VM are similar using first index 2, 3.
!
!    Output, integer ( kind = 4 ) NTRIF(1:NPOLH), the number of triangular faces in
!    triangulation of each polyhedron.
!
!    Output, integer ( kind = 4 ) NTETRA(1:NPOLH), the number of tetrahedra in
!    triangulation of each polyhedron.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Worksapce, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxbt
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) maxht
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxvm
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npolh
  integer ( kind = 4 ) nvert

  integer ( kind = 4 ) btst(nface+1)
  integer ( kind = 4 ) btl(3,maxbt)
  integer ( kind = 4 ) ceang
  integer ( kind = 4 ) cfvl
  integer ( kind = 4 ) chvl
  integer ( kind = 4 ) cnrml
  integer ( kind = 4 ) crit
  real    ( kind = 8 ) ctrx
  real    ( kind = 8 ) ctry
  real    ( kind = 8 ) ctrz
  real    ( kind = 8 ) diamsq
  integer ( kind = 4 ) drem
  real    ( kind = 8 ) eang(nvert)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edno(nvert)
  integer ( kind = 4 ) edst(nvert)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,nface)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fcst(nface+1)
  integer ( kind = 4 ) fvl(6,nvert)
  integer ( kind = 4 ) g
  real    ( kind = 8 ) h(npolh)
  integer ( kind = 4 ) hfl(npolh)
  integer ( kind = 4 ) ht(maxht)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifc
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) ivm
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jv
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kfc
  integer ( kind = 4 ) kht
  integer ( kind = 4 ) kvm
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lk
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) maxsv
  integer ( kind = 4 ) n
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbmv
  integer ( kind = 4 ) ncface
  integer ( kind = 4 ) ncvert
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) nimv
  integer ( kind = 4 ) npt
  real    ( kind = 8 ) nrml(3,nface)
  integer ( kind = 4 ) nsface
  integer ( kind = 4 ) nsvc
  integer ( kind = 4 ) nsvert
  integer ( kind = 4 ) ntetra(npolh)
  integer ( kind = 4 ) ntrif(npolh)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,*)
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) prime
  logical              rflag
  integer ( kind = 4 ) sf
  integer ( kind = 4 ) sfvl
  real    ( kind = 8 ) shrf
  integer ( kind = 4 ) shvl
  integer ( kind = 4 ) sizht
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) svcl
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
  integer ( kind = 4 ) vm(maxvm)
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xfhv(3,npolh+1)
!
!  Determine NCFACE = maximum number of faces in a polyhedron and NCVERT = 2 *
!  maximum number of edges in a polyhedron.  Array FCST is temporarily used.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  ncface = 0
  ncvert = 0

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

    fcst(f) = k
  end do

  do p = 1,npolh

    j = 0
    k = 0
    i = hfl(p)

    do

      f = abs(pfl(1,i))
      j = j + 1
      k = k + fcst(f)
      i = pfl(2,i)

      if ( i == hfl(p) ) then
        exit
      end if

    end do

    ncface = max ( ncface, j )
    ncvert = max ( ncvert, k )

  end do
!
!  Generate 2D Delaunay triangulations on convex boundary faces
!  of polyhedra.  Then generate interior mesh vertices and tetrahedra
!  for each convex polyhedron.  Start of IWK is used for mapping
!  global vertex indices to local indices for a polyhedron.
!
  call tribfc(h,nvc,nface,nvert,maxvc,maxbt,maxiw,maxwk,vcl,facep, &
     nrml,fvl,edst,edno,fcst,btst,btl,iwk,wk, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  chvl = nvc + 1
  cfvl = chvl + ncface
  irem = cfvl + 5*ncvert
  cnrml = 1
  ceang = cnrml + 3*ncface
  drem = ceang + ncvert

  if ( maxiw < irem - 1 ) then
    ierr = 6
    return
  else if ( maxwk < drem - 1 ) then
    ierr = 7
    return
  end if

  iwk(1:nvc) = 0
  kfc = 0
  kht = 0
  kvm = 0
  xfhv(1,1) = 1
  xfhv(2,1) = 1
  xfhv(3,1) = 1
  rflag = .true.

  do p = 1, npolh

    call dsconv(p,hfl(p),facep,nrml,fvl,eang,pfl,ncface,ncvert, &
      iwk(chvl),wk(cnrml),iwk(cfvl),wk(ceang))

    irem = cfvl + 5*ncvert

    if ( maxiw < irem + ncvert - 1 ) then
      ierr = 6
      return
    end if

    call rmcpfc(ncface,ncvert,iwk(chvl),wk(cnrml),iwk(cfvl), &
      wk(ceang),iwk(irem))

    call rmcled(ncface,ncvert,iwk(chvl),iwk(cfvl))

    maxsv = 2*ncvert
    shvl = cfvl + 5*ncvert
    sfvl = shvl + ncface
    irem = sfvl + 5*maxsv
    svcl = ceang + ncvert
    drem = svcl + 3*maxsv

    if ( maxiw < irem ) then
      ierr = 6
      return
    else if ( maxwk < drem ) then
      ierr = 7
      return
    end if

    call shrnk3(shrf*h(p),ncface,vcl,iwk(chvl),wk(cnrml),iwk(cfvl), &
      wk(ceang),maxsv,maxiw-irem+1,maxwk-drem+1,nsvc,nsface, &
      nsvert,wk(svcl),iwk(shvl),iwk(sfvl),iwk(irem),wk(drem), ierr )

    if ( ierr == 13) then
      ierr = 330
    end if

    if ( ierr /= 0 ) then
      return
    end if

    if ( 4 <= nsface ) then
      call diam3(nsvc,wk(svcl),ibot,itop,diamsq)
      drem = svcl + 3*nsvc
      call intmvg(nsvc,nsface,nsvert,wk(svcl),iwk(shvl),iwk(sfvl), &
        ibot,itop,h(p),maxvc-nvc,maxwk-drem+1,nimv,vcl(1,nvc+1), &
        wk(drem), ierr )

      if ( ierr /= 0 ) then
        return
      end if

      nvc = nvc + nimv
    else
      nimv = 0
    end if

    nbmv = 0
    bf_num = 0
!
!  Add vertices of convex polyhedron to VM.
!
    ivm = kvm
    i = hfl(p)

60  continue

    f = abs(pfl(1,i))
    j = facep(1,f)

!70  continue

    do

      lj = fvl(loc,j)

      if ( iwk(lj) == 0 ) then

        ivm = ivm + 1

        if ( maxvm < ivm ) then
          ierr = 20
          return
        end if

        vm(ivm) = lj
        iwk(lj) = ivm - kvm
        nbmv = nbmv + 1

      end if

      j = fvl(succ,j)

      if ( j == facep(1,f) ) then
        exit
      end if

    end do

    i = pfl(2,i)

    if ( i /= hfl(p)) go to 60

    lv = ivm
!
!  Add mesh vertices on interior of edges and faces to VM.
!
    i = hfl(p)

80  continue

    sf = pfl(1,i)
    f = abs(sf)
    bf_num = bf_num + (btst(f+1) - btst(f))
    j = facep(1,f)
    lj = fvl(loc,j)

90  continue

    k = fvl(succ,j)
    lk = fvl(loc,k)

    if ( 0 < (lk - lj)*sf ) then
      g = fvl(facn,fvl(edgc,j))
    else
      g = fvl(facn,fvl(edga,j))
    end if

    if ( f < g .and. 0 < edst(j) ) then

      n = abs(edno(j))

      if ( maxvm < ivm + n ) then
        ierr = 20
        return
      end if

      nbmv = nbmv + n

      do l = edst(j),edst(j)+n-1
        ivm = ivm + 1
        vm(ivm) = l
        iwk(l) = ivm - kvm
      end do

    end if

    j = k
    lj = lk

    if ( j /= facep(1,f)) go to 90

    n = fcst(f+1) - fcst(f)

    if ( 0 < n ) then

      if ( maxvm < ivm + n ) then
        ierr = 20
        return
      end if

      nbmv = nbmv + n

      do l = fcst(f),fcst(f)+n-1
        ivm = ivm + 1
        vm(ivm) = l
        iwk(l) = ivm - kvm
      end do

    end if

    i = pfl(2,i)

    if ( i /= hfl(p) ) go to 80
!
!  Add mesh vertices in interior of polyhedron to VM.
!
    if ( 0 < nimv ) then

      if ( maxvm < ivm + nimv ) then
        ierr = 20
        return
      end if

      do l = nvc-nimv+1,nvc
        ivm = ivm + 1
        vm(ivm) = l
      end do

    else

      ivm = ivm + 1

      if ( maxvm < ivm ) then
        ierr = 20
        return
      else if ( maxvc <= nvc ) then
        ierr = 14
        return
      end if

      vm(ivm) = nvc + 1
      ctrx = 0.0D+00
      ctry = 0.0D+00
      ctrz = 0.0D+00
      jv = kvm + 1

      do i = jv, ivm-1
        j = vm(i)
        ctrx = ctrx + vcl(1,j)
        ctry = ctry + vcl(2,j)
        ctrz = ctrz + vcl(3,j)
      end do

      ctrx = ctrx / real ( nbmv, kind = 8 )
      ctry = ctry / real ( nbmv, kind = 8 )
      ctrz = ctrz / real ( nbmv, kind = 8 )

      j = vm(jv)

      vcl(1,nvc+1) = 0.5D+00 * ( ctrx + vcl(1,j) )
      vcl(2,nvc+1) = 0.5D+00 * ( ctry + vcl(2,j) )
      vcl(3,nvc+1) = 0.5D+00 * ( ctrz + vcl(3,j) )

      rflag = .true.

    end if

    npt = ivm - kvm
    sizht = prime(3*npt/2)

    if ( sizht < 3*npt/2 ) then
      ierr = 100
      return
    else if ( maxht < kht + sizht ) then
      ierr = 21
      return
    else if ( fc_max < kfc + bf_num ) then
      ierr = 11
      return
    end if
!
!  Set up boundary faces with local vertex indices.
!
    ifc = kfc
    i = hfl(p)

!140 continue

    do

      f = abs(pfl(1,i))

      do k = btst(f), btst(f+1)-1
        ifc = ifc + 1
        fc(1:3,ifc) = iwk(btl(1:3,k))
      end do

      i = pfl(2,i)

      if ( i == hfl(p) ) then
        exit
      end if

    end do

    do i = kvm+1, kvm+nbmv
      iwk(vm(i)) = 0
    end do
!
!  Construct boundary-constrained triangulations. Dummy array IWK
!  is used for BF in call to IMPTR3, since BF is not referenced.
!
180 continue

    call bcdtri(rflag,bf_num,nbmv,nimv,sizht,fc_max-kfc,vcl,vm(kvm+1), &
      nfc,ntrif(p),ntetra(p),fc(1,kfc+1),ht(kht+1), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( nimv == 0 .and. 0 < vm(ivm) ) then

      if ( jv < lv ) then
        jv = jv + 1
        j = vm(jv)
        vcl(1,nvc+1) = 0.5D+00 * ( ctrx + vcl(1,j) )
        vcl(2,nvc+1) = 0.5D+00 * ( ctry + vcl(2,j) )
        vcl(3,nvc+1) = 0.5D+00 * ( ctrz + vcl(3,j) )
        go to 180
      else if ( jv == lv ) then
        jv = jv + 1
        vcl(1,nvc+1) = ctrx
        vcl(2,nvc+1) = ctry
        vcl(3,nvc+1) = ctrz
        rflag = .false.
        go to 180
      else
        nvc = nvc + 1
      end if

    end if

    if ( 1 <= crit .and. crit <= 3 ) then

      call imptr3 ( .true., .true., &
        crit, npt, sizht, fc_max-kfc, vcl, vm(kvm+1), nfc, ntetra(p), iwk, &
        fc(1,kfc+1), ht(kht+1), ntrif(p), ierr )

    end if

    if ( ierr /= 0 ) then
      return
    end if

    kfc = kfc + nfc
    kht = kht + sizht
    kvm = ivm
    xfhv(1,p+1) = kfc + 1
    xfhv(2,p+1) = kht + 1
    xfhv(3,p+1) = kvm + 1

  end do

  return
end
