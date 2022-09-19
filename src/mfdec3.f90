subroutine mfdec3 ( hflag, umdf, kappa, angacc, angedg, dmin, nmin, ntetd, &
  nsflag, nvc, nface, nvert, npolh, npf, maxvc, maxfp, maxfv, maxhf, maxpf, &
  maxiw, maxwk, vcl, facep, factyp, nrml, fvl, eang, hfl, pfl, ivrt, xivrt, &
  ifac, xifac, wid, facval, edgval, vrtval, vol, psi, htsiz, maxedg, ht, &
  edge, listev, infoev, iwk, wk, ierr )

!*****************************************************************************80
!
!! MFDEC3 subdivides polyhedra to control the mesh distribution function.
!
!  Discussion:
!
!    This routine further subdivides convex polyhedra so that the variation
!    of heuristic or user-supplied mesh distribution function in
!    each polyhedron is limited.
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
!    Input, external real ( kind = 8 ) UMDF(X,Y,Z), the user-supplied mdf
!    with d.p arguments.
!
!    Input, real ( kind = 8 ) KAPPA, the mesh smoothness parameter in
!    interval [0.0,1.0], used iff HFLAG is TRUE.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by cut faces.
!
!    Input, real ( kind = 8 ) ANGEDG, the angle parameter in radians used
!    to determine allowable points on edges as possible endpoints of edges
!    of cut faces.
!
!    Input, real ( kind = 8 ) DMIN, the parameter used to determine if
!    variation of mdf in
!    Input, NMIN, the parameter used to determine if 'sufficiently large'
!    number of tetrahedra in polyhedron.
!
!    Input, integer ( kind = 4 ) NTETD, the desired number of tetrahedra in mesh.
!
!    Input, logical NSFLAG, TRUE if continue to next polyhedron when no
!    separator face is found for a polyhedron, FALSE if terminate with
!    error 336 when no separator face is found for a polyhedron.
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
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    greater than or equal to 2*(NE + MAX ( NF, NV ) ) where NE, NF, NV
!    are maximum number of edges, faces, vertices in any polyhedron of updated
!    decomposition.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    greater than or equal to MAX ( NPOLH, NE+MAX ( 2*NF, 3*NV ) ) where NPOLH
!    is input value.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list: row 1
!    is head pointer, rows 2 and 3 are signed polyhedron indices.
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
!    Input/output, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices
!    in PFL for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices
!    for each polyhedron; row 2 used for link.
!
!    Input, integer ( kind = 4 ) IVRT(1:NVERT), XIVRT(1:NFACE+1), IFAC(1:NPF),
!    XIFAC(1:NPOLH+1), real ( kind = 8 ) WID(1:NPOLH), FACVAL(1:NFACE),
!    EDGVAL(1:NVERT), VRTVAL(1:NVC), arrays output from routine DSMDF3 if
!    HFLAG is .TRUE.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime
!    number which is greater than or equal to the number of vertices in
!    a polyhedron of decomposition.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array;
!    should be greater than or equal to the maximum number of edges in a
!    polyhedron of decomposition.
!
!    Ouptut, real ( kind = 8 ) VOL(1:NPOLH), the volume of convex polyhedra
!    in decomposition.
!
!    Output, real ( kind = 8 ) PSI(1:NPOLH), the mean mdf values in the
!    convex polyhedra.
!
!    Workspace, integer HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG), the hash table and
!    edge records used to determine entries of LISTEV.
!
!    Workspace, integer LISTEV(1:*), used by routines PRMDF3, INTPH; size
!    must be greater than or equal to NCFACE+NCEDGE+NCVC where NCFACE =
!    maximum number of faces in a polyhedron (of input decomposition),
!    NCEDGE = maximum number of edges in a polyhedron,
!    NCVC = maximum number of vertices in a polyhedron.
!
!    Workspace, integer INFOEV(1:4,1:*), used by routines PRMDF3, INTPH;
!    size must be greater than or equal to NCFACE+NCEDGE.
!
!    [Note: It is assumed there is enough space for the arrays.
!    HT, EDGE, LISTEV and INFOEV are needed only if HFLAG is .TRUE.]
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxhf
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert

  logical              aflag
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) angedg
  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) c1
  real    ( kind = 8 ) c2
  integer ( kind = 4 ) ccw
  integer ( kind = 4 ) cdang
  integer ( kind = 4 ) cedge
  real    ( kind = 8 ) dmin
  real    ( kind = 8 ) dtol
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edge(4,maxedg)
  real    ( kind = 8 ) edgval(nvert)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  real    ( kind = 8 ) facval(nface)
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) hfl(maxhf)
  logical              hflag
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac(npf)
  integer ( kind = 4 ) indf
  real    ( kind = 8 ) infoev(4,*)
  real    ( kind = 8 ) intreg
  integer ( kind = 4 ) ivrt(nvert)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) kappa
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) l
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) listev(*)
  integer ( kind = 4 ) ll
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lp
  real    ( kind = 8 ) mdfint
  real    ( kind = 8 ) mean
  integer ( kind = 4 ) meanf
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mxcos
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nce
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nedev
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nfcev
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) np
  real    ( kind = 8 ) nrml(3,maxfp)
  real    ( kind = 8 ) nrmlc(4)
  logical              nsflag
  integer ( kind = 4 ) ntetd
  integer ( kind = 4 ) nvcin
  integer ( kind = 4 ) nvrev
  integer ( kind = 4 ) pfl(2,maxpf)
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) psi(maxhf)
  real    ( kind = 8 ) stdv
  integer ( kind = 4 ) stdvf
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), external :: umdf
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) vol(maxhf)
  real    ( kind = 8 ) volreg
  real    ( kind = 8 ) vrtval(nvc)
  real    ( kind = 8 ) wid(npolh)
  real    ( kind = 8 ) widp
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xifac(npolh+1)
  integer ( kind = 4 ) xivrt(nface+1)
!
!  WK(1:NPOLH) is used for mdf standard deviation in polyhedra.
!  Compute VOLREG = volume of region and INTREG = estimated integral
!  of MDF3(X,Y) or UMDF(X,Y).
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( maxwk < npolh ) then
    ierr = 7
    return
  end if

  nvcin = nvc
  volreg = 0.0D+00
  intreg = 0.0D+00
  widp = 0.0D+00

  do i = 1, npolh

    if ( hflag ) then
      widp = wid(i)
      call prmdf3(i,widp,nvcin,vcl,nrml,ivrt,xivrt,ifac,xifac, &
        facval,edgval,vrtval,nfcev,nedev,nvrev,listev,infoev, &
        htsiz,maxedg,ht,edge, ierr )
      if ( ierr /= 0 ) then
        return
      end if
    end if

    n = 0
    cntr(1:3) = 0.0D+00
    j = hfl(i)

10  continue

    f = abs ( pfl(1,j) )
    k = facep(1,f)

    do

      l = fvl(loc,k)
      cntr(1:3) = cntr(1:3) + vcl(1:3,l)
      n = n + 1
      k = fvl(succ,k)

      if ( k == facep(1,f) ) then
        exit
      end if

    end do

    j = pfl(2,j)
    if ( j /= hfl(i)) then
      go to 10
    end if

    cntr(1:3) = cntr(1:3) / real ( n, kind = 8 )

    call intph(hflag,umdf,hfl(i),widp,nfcev,nedev,nvrev,listev, &
      infoev,ivrt,facval,edgval,vrtval,vcl,facep,fvl,pfl,cntr, &
      mdfint,psi(i),wk(i),vol(i),1,iwk,wk,wk)
    volreg = volreg + vol(i)
    intreg = intreg + mdfint

  end do
!
!  If HFLAG, compute mean mdf values from KAPPA, etc. Scale PSI(I)'s
!  so that integral in region is 1. Determine which polyhedra need to
!  be further subdivided (indicated by negative PSI(I) value).
!
  if ( hflag ) then
    c1 = (1.0D+00 - kappa)/intreg
    c2 = kappa/volreg
  else
    c1 = 1.0D+00/intreg
    c2 = 0.0D+00
  end if

  do i = 1, npolh
    psi(i) = psi(i)*c1 + c2
    if ( psi(i)*dmin < c1 * wk(i) ) then
      if ( nmin < ntetd*psi(i)*vol(i) ) then
         psi(i) = -psi(i)
      end if
    end if
  end do
!
!  Further subdivide polygons for which DMIN < STDV/MEAN and
!  (estimated number of tetrahedra) greater than NMIN.
!
  mxcos = cos(angedg)
  np = npolh
  cedge = 1
  cdang = 1

  do i = 1, np

    if ( 0.0D+00 <= psi(i) ) then
      cycle
    end if

    if ( hflag ) then
      widp = wid(i)
      call prmdf3(i,widp,nvcin,vcl,nrml,ivrt,xivrt,ifac,xifac, &
        facval,edgval,vrtval,nfcev,nedev,nvrev,listev,infoev, &
        htsiz,maxedg,ht,edge, ierr )
      if ( ierr /= 0 ) then
        return
      end if
    end if

    lp = npolh + 1
    kp = i

50  continue

    nf = 0
    n = 0
    cntr(1:3) = 0.0D+00
    j = hfl(kp)

60  continue

    nf = nf + 1
    f = abs ( pfl(1,j) )
    k = facep(1,f)

70  continue

    l = fvl(loc,k)
    cntr(1:3) = cntr(1:3) + vcl(1:3,l)
    n = n + 1
    k = fvl(succ,k)
    if ( k /= facep(1,f)) then
      go to 70
    end if

    j = pfl(2,j)

    if ( j /= hfl(kp)) then
      go to 60
    end if

    cntr(1:3) = cntr(1:3)/ real ( n, kind = 8 )
    ne = n/2
    indf = cedge + n
    meanf = ne
    stdvf = meanf + nf

    if ( maxiw < indf + nf + nf - 1 ) then
      ierr = 6
      return
    else if ( maxwk < stdvf + nf - 1 ) then
      ierr = 7
      return
    end if

    call intph(hflag,umdf,hfl(kp),widp,nfcev,nedev,nvrev, &
      listev,infoev,ivrt,facval,edgval,vrtval,vcl,facep,fvl, &
      pfl,cntr,mdfint,mean,stdv,vol(kp),nf,iwk(indf), &
      wk(meanf),wk(stdvf))

    psi(kp) = mean*c1 + c2

    if ( psi(kp) * dmin < c1 * stdv ) then

      if ( nmin < ntetd * psi(kp) * vol(kp) ) then

        sum2 = 0.0D+00
        j = hfl(kp)

80      continue

        f = pfl(1,j)

        if ( 0 < f ) then
          ccw = succ
        else
          ccw = pred
          f = -f
        end if

        k = facep(1,f)
        l = fvl(loc,k)

90      continue

        kk = fvl(ccw,k)
        ll = fvl(loc,kk)

        if ( l < ll ) then
          leng = sqrt((vcl(1,l)-vcl(1,ll))**2+(vcl(2,l) &
            -vcl(2,ll))**2 + (vcl(3,l)-vcl(3,ll))**2)
          sum2 = sum2 + leng
        end if

        k = kk
        l = ll

        if ( k /= facep(1,f)) then
          go to 90
        end if

        j = pfl(2,j)
        if ( j /= hfl(kp)) then
          go to 80
        end if

        dtol = tol * sum2 / real ( ne, kind = 8 )

        call sfc1mf(kp,cntr,mean,nf,iwk(indf),wk(meanf), &
          angacc,mxcos,dtol,nvc,maxvc,vcl,facep,nrml,fvl, &
          eang,nrmlc,nce,iwk(cedge),wk(cdang),aflag, &
          iwk(indf+nf), ierr )

        if ( ierr /= 0 ) then
          return
        end if

        if ( aflag ) then
          go to 110
        end if

        stdv = -1.0D+00

        do j = 0, nf-1

          t = wk(stdvf+j) / wk(meanf+j)

          if ( stdv < t ) then
            stdv = t
            k = j
          end if

        end do

        f = iwk(indf+k)
        n = ne + 2 - nf

        if ( maxiw < indf + n + n - 1 ) then
          ierr = 6
          return
        else if ( maxwk < meanf + 3*n - 1 ) then
          ierr = 7
          return
        end if

        call sfc2mf(kp,f,hflag,umdf,widp,nfcev,nedev,nvrev, &
          listev,infoev,ivrt,facval,edgval,vrtval,cntr, &
          angacc,angedg,mxcos,dtol,nvc,maxvc,vcl,facep,nrml, &
          fvl,eang,nrmlc,nce,iwk(cedge),wk(cdang),aflag, &
          iwk(indf),wk(meanf), ierr )

        if ( ierr /= 0 ) then
          return
        end if

        if ( aflag ) then
          go to 110
        end if

        call sfcshp(kp,hfl(kp),cntr,angacc,mxcos,dtol,nvc, &
          maxvc,vcl,facep,nrml,fvl,eang,pfl,nrmlc,nce, &
          iwk(cedge),wk(cdang),aflag,n,iwk(indf),wk(meanf), ierr )

        if ( ierr /= 0 ) then
          return
        end if

        if ( .not. aflag ) then
          if ( nsflag ) then
            if ( msglvl == 4 ) then
              write ( *,600)
            end if
            go to 120
          end if
          ierr = 336
          return
        end if

110     continue

        call insfac(kp,nrmlc,nce,iwk(cedge),wk(cdang),nvc, &
          nface,nvert,npolh,npf,maxfp,maxfv,maxhf,maxpf,vcl, &
          facep,factyp,nrml,fvl,eang,hfl,pfl,ierr)

        if ( ierr /= 0 ) then
          return
        end if

        psi(kp) = -psi(kp)
        psi(npolh) = psi(kp)

      end if

    end if

    if ( psi(kp) < 0.0D+00 ) then
      go to 50
    end if

120 continue

    if ( kp == i ) then
      kp = lp
    else
      kp = kp + 1
    end if

    if ( kp <= npolh ) then
      go to 50
    end if

  end do

  600 format (4x,'*** no separator face found')

  return
end
