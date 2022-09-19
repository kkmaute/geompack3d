subroutine eqdis3 ( hflag, umdf, kappa, angacc, angedg, dmin, nmin, ntetd, &
  nsflag, nvc, nface, nvert, npolh, npf, maxvc, maxfp, maxfv, maxhf, maxpf, &
  maxiw, maxwk, vcl, facep, factyp, nrml, fvl, eang, hfl, pfl, vol, psi, h, &
  iwk, wk, ierr )

!*****************************************************************************80
!
!! EQDIS3 subdivides polyhedra for equidistribution.
!
!  Discussion:
!
!    This routine further subdivides a convex polyhedra so that an approximate
!    equidistributing tetrahedral mesh can be constructed with
!    respect to heuristic or user-supplied mesh distribution
!    function, and determine tetrahedron size for each polyhedron
!    of decomposition.
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
!    Input, external real ( kind = 8 ) UMDF(X,Y,Z), d.p user-supplied mdf with
!    d.p arguments.
!
!    Input, real ( kind = 8 ) KAPPA, the mesh smoothness parameter in
!    interval [0.0,1.0], used iff HFLAG is TRUE.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by cut faces.
!
!    Input, real ( kind = 8 ) ANGEDG, the angle parameter in radians used to
!    determine allowable points on edges as possible endpoints of edges of
!    cut faces.
!
!    Input, real ( kind = 8 ) DMIN, the parameter used to determine if
!    variation of mdf in polyhedron is 'sufficiently high'.
!
!    Input, integer ( kind = 4 ) NMIN, the parameter used to determine if 'sufficiently
!    large' number of tetrahedra in polyhedron.
!
!    Input, integer ( kind = 4 ) NTETD, the desired number of tetrahedra in mesh.
!
!    Input, logical NSFLAG, the TRUE if continue to next polyhedron when no
!    separator face is found for a polyhedron, FALSE if terminate with
!    error 336 when no separator face is found for a polyhedron.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL.
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces or positions used
!    in FACEP array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input/output, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions
!    used in HFL array.
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
!    Input, integer ( kind = 4 ) MAXHF, the maximum size available for HFL array.
!
!    Input, integer ( kind = 4 ) MAXPF, the maximum size available for PFL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array;
!    should be greater than or equal to
!    MAX ( S*(NVERT+NFACE+NPF+NPOLH+2 + 2*NCVC+5*NCEDGE+NCFACE)
!    + 2*(NE + MAX ( NF, NV ) ), 2*NPOLHout ) where
!    S = 1 or 0 if HFLAG is TRUE or FALSE, NVERT to NPOLH are
!    input values, NCVC = max no. of vertices in a polyhedron (of
!    input decomposition), NCEDGE = max no. of edges in a
!    polyhedron, NCFACE = max number of faces in a polyh, NE,NF,NV
!    are max number of edges, faces, vertices in any polyhedron
!    of updated decomposition, NPOLHout is output value of NPOLH.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    greater than or equal to
!    S*(NPOLH+NFACE+NVERT+NVC + 4*NCFACE+4*NCEDGE) +
!    MAX ( NPOLH, NE+MAX ( 2*NF, 3*NV ) ) where NVC is input value.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate
!    list.
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
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list;
!    see routine DSPHDC.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the angles at edges
!    common to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J),
!    determined by EDGC field.
!
!    Input/output, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices
!    in PFL for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face
!    indices for each polyhedron; row 2 used for link.
!
!    Output, real ( kind = 8 ) VOL(1:NPOLH), the volume of convex polyhedra
!    in decomposition.
!
!    Output, real ( kind = 8 ) PSI(1:NPOLH), the mean mdf values in the
!    convex polyhedra.
!
!    Output, real ( kind = 8 ) H(1:NPOLH), the tetrahedron size for
!    convex polyhedra.
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
  integer ( kind = 4 ) maxhf
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) angedg
  real    ( kind = 8 ) dmin
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edgval
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) facval
  integer ( kind = 4 ) fvl(6,maxfv)
  real    ( kind = 8 ) h(maxhf)
  integer ( kind = 4 ) hfl(maxhf)
  logical              hflag
  integer ( kind = 4 ) ht
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) infoev
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(maxiw)
  real    ( kind = 8 ) kappa
  integer ( kind = 4 ) listev
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncedge
  integer ( kind = 4 ) ncface
  integer ( kind = 4 ) ncvc
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  real    ( kind = 8 ) nrml(3,maxfp)
  logical              nsflag
  integer ( kind = 4 ) ntetd
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) pfl(2,maxpf)
  integer ( kind = 4 ) prime
  real    ( kind = 8 ) psi(maxhf)
  real    ( kind = 8 ), external :: umdf
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) vol(maxhf)
  integer ( kind = 4 ) vrtval
  integer ( kind = 4 ) wid
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xifac
  integer ( kind = 4 ) xivrt

  ierr = 0
  ivrt = 1
  wid = 1

  if ( hflag ) then

    xivrt = ivrt + nvert
    ifac = xivrt + nface + 1
    xifac = ifac + npf
    m = xifac + npolh
    facval = wid + npolh
    edgval = facval + nface
    vrtval = edgval + nvert
    n = vrtval + nvc - 1

    if ( maxiw < m ) then
      ierr = 6
      return
    else if ( maxwk < n ) then
      ierr = 7
      return
    end if

    call dsmdf3(nvc,nface,nvert,npolh,maxiw-m,maxwk-n,vcl,facep, &
      nrml,fvl,eang,hfl,pfl,iwk(ivrt),iwk(xivrt),iwk(ifac), &
      iwk(xifac),wk(wid),wk(facval),wk(edgval),wk(vrtval),ncface, &
      ncedge,iwk(m+1),wk(n+1),ierr)

    if ( ierr /= 0 ) then
      return
    end if

    ncvc = ncedge - 2
    htsiz = prime(ncvc)
    ht = m + 1
    edge = ht + htsiz
    listev = edge + 4 * ncedge
    m = listev + ncface + ncedge + ncvc - 1
    infoev = n + 1

    n = infoev + 4 * ( ncface + ncedge ) - 1

    if ( maxiw < m ) then
      ierr = 6
      return
    else if ( maxwk < n ) then
      ierr = 7
      return
    end if

  else

    xivrt = 1
    ifac = 1
    xifac = 1
    facval = 1
    edgval = 1
    vrtval = 1
    htsiz = 1
    ncedge = 1
    ht = 1
    edge = 1
    listev = 1
    infoev = 1
    m = 0
    n = 0

  end if

  call mfdec3(hflag,umdf,kappa,angacc,angedg,dmin,nmin,ntetd,nsflag, &
    nvc,nface,nvert,npolh,npf,maxvc,maxfp,maxfv,maxhf,maxpf, &
    maxiw-m,maxwk-n,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
    iwk(ivrt),iwk(xivrt),iwk(ifac),iwk(xifac),wk(wid),wk(facval), &
    wk(edgval),wk(vrtval),vol,psi,htsiz,ncedge,iwk(ht),iwk(edge), &
    iwk(listev),wk(infoev),iwk(m+1),wk(n+1), ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( maxiw < npolh + npolh ) then
    ierr = 6
    return
  end if

  call tetsiz ( ntetd, npolh, facep, hfl, pfl, vol, psi, h, iwk, iwk(npolh+1) )

  return
end
