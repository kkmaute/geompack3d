subroutine sfc2mf ( p, f, hflag, umdf, widp, nfcev, nedev, nvrev, listev, &
  infoev, ivrt, facval, edgval, vrtval, cntr, angacc, angedg, mxcos, dtol, &
  nvc, maxvc, vcl, facep, nrml, fvl, eang, nrmlc, nce, cedge, cdang, aflag, &
  indv, val, ierr )

!*****************************************************************************80
!
!! SFC2MF finds a separator or cut face in a convex polyhedron.
!
!  Discussion:
!
!    This routine attempts to find a separator or cut face in a convex
!    polyhedron P based on a mesh distribution function by starting with an edge
!    in the interior of face F which separates higher mdf values
!    from lower ones on the face.
!
!    Accept a cut face if it creates no small angles or short edges.
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
!    Input, integer ( kind = 4 ) P, the polyhedron index.
!
!    Input, integer ( kind = 4 ) F, the face index.
!
!    Input, logical HFLAG, is TRUE if heuristic MDF, FALSE if user-supplied
!    MDF.
!
!    Input, external real ( kind = 8 ) UMDF(X,Y,Z), the user-supplied mdf
!    with d.p arguments.
!
!    Input, real ( kind = 8 ) WIDP, the width of original polyhedron
!    of decomposition.
!
!    Input, integer ( kind = 4 ) NFCEV, NEDEV, NVREV, LISTEV(1:NFCEV+NEDEV+NVREV),
!    INFOEV(1:4,1:NFCEV+NEDEV), output from routine PRMDF3.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), real ( kind = 8 ) FACVAL(1:*),
!    EDGVAL(1:*), VRTVAL(1:*), arrays output from routine DSMDF3.
!
!    [Note: Parameters WIDP to VRTVAL are used only if HFLAG = TRUE]
!
!    Input, real ( kind = 8 ) CNTR(1:3), the weighted centroid of polyhedron.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral
!    angle in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) ANGEDG, the angle parameter in radians used
!    to determine allowable points on edges as possible endpoints of edges
!    of cut faces.
!
!    Input, real ( kind = 8 ) MXCOS, the maximum cosine of angle allowed
!    for angles subtended by new subedges with respect to
!    centroid = COS(ANGEDG); exception is that endpoints of start edge may
!    be midpoints of existing edges.
!
!    Input, real ( kind = 8 ) DTOL, the absolute tolerance to determine
!    if point is on cut plane.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC).  On input, vertex
!    coordinate list.  On output, some temporary or permanent entries
!    may be added at end of array.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:*), the edge angles.
!
!    Output, real ( kind = 8 ) NRMLC(1:4), the unit normal vector of cut
!    plane plus right hand side constant term of plane equation (if
!    acceptable).
!
!    Output, integer ( kind = 4 ) NCE, the number of edges in cut face (if acceptable); it
!    is assumed there is enough space in the following two arrays.
!
!    Output, integer ( kind = 4 ) CEDGE(1:2,0:NCE), real ( kind = 8 ) CDANG(1:NCE),
!    information describing cut polygon as output by routine SEPFAC.
!
!    Output, logical AFLAG, TRUE iff separator face is found and acceptable.
!
!    Workspace, integer INDV(1:2*NE), where NE is number of edges on face F.
!
!    Workspace, real ( kind = 8 ) VAL(1:2*NE).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxvc

  integer ( kind = 4 ) a
  logical              aflag
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) angedg
  integer ( kind = 4 ) b
  real    ( kind = 8 ) cdang(*)
  real    ( kind = 8 ) ce
  integer ( kind = 4 ) cedge(2,0:*)
  real    ( kind = 8 ) cn
  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) cosed2
  real    ( kind = 8 ) cosmax
  real    ( kind = 8 ) cosmin
  real    ( kind = 8 ) da
  real    ( kind = 8 ) db
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dtol
  real    ( kind = 8 ) eang(*)
  real    ( kind = 8 ) edg(3)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) edgval(*)
  real    ( kind = 8 ) enr(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  real    ( kind = 8 ) facval(*)
  real    ( kind = 8 ), parameter, dimension ( 3 ) :: fract= (/ &
    0.5D+00, 0.4D+00, 0.6D+00 /)
  integer ( kind = 4 ) fvl(6,*)
  logical              hflag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii(2)
  integer ( kind = 4 ) indv(*)
  real    ( kind = 8 ) infoev(4,*)
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj(2)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kvc
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) listev(*)
  integer ( kind = 4 ) ll(2)
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mdf3
  integer ( kind = 4 ) mm(2)
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mxcos
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nce
  integer ( kind = 4 ) nedev
  logical              neg
  integer ( kind = 4 ) nfcev
  real    ( kind = 8 ) nrml(3,*)
  real    ( kind = 8 ) nrmlc(4)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvrev
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) s(2)
  integer ( kind = 4 ) sf
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) umdf
  real    ( kind = 8 ) va(3)
  real    ( kind = 8 ) val(*)
  real    ( kind = 8 ) vave
  real    ( kind = 8 ) vb(3)
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) widp
!
!  Evaluate MDF at vertices of face and midpoints of edges
!  subtending large angles with respect to centroid.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  aflag = .false.
  cosed2 = cos ( 2.0D+00 * angedg )
  kvc = nvc
  n = 0
  a = facep(1,f)
  la = fvl(loc,a)
  va(1) = vcl(1,la) - cntr(1)
  va(2) = vcl(2,la) - cntr(2)
  va(3) = vcl(3,la) - cntr(3)
  da = va(1)**2 + va(2)**2 + va(3)**2

10 continue

  b = fvl(succ,a)
  lb = fvl(loc,b)
  vb(1) = vcl(1,lb) - cntr(1)
  vb(2) = vcl(2,lb) - cntr(2)
  vb(3) = vcl(3,lb) - cntr(3)
  db = vb(1)**2 + vb(2)**2 + vb(3)**2
  dotp = (va(1)*vb(1) + va(2)*vb(2) + va(3)*vb(3))/sqrt(da*db)

  if ( dotp < cosed2 ) then

    kvc = kvc + 1

    if ( maxvc < kvc ) then
      ierr = 14
      return
    end if

    vcl(1:3,kvc) = 0.5D+00*(vcl(1:3,la) + vcl(1:3,lb))

    n = n + 1
    indv(n) = -kvc

    if ( hflag ) then
      val(n) = mdf3(vcl(1,kvc),vcl(2,kvc),vcl(3,kvc),widp, &
        nfcev,nedev,nvrev,listev,infoev,ivrt,facval,edgval, &
        vrtval,vcl)
    else
      val(n) = umdf(vcl(1,kvc),vcl(2,kvc),vcl(3,kvc))
    end if

  end if

  n = n + 1
  indv(n) = b

  if ( hflag ) then
    val(n) = mdf3(vcl(1,lb),vcl(2,lb),vcl(3,lb),widp,nfcev, &
      nedev,nvrev,listev,infoev,ivrt,facval,edgval,vrtval,vcl)
  else
    val(n) = umdf(vcl(1,lb),vcl(2,lb),vcl(3,lb))
  end if

  a = b

  if ( a /= facep(1,f) ) then
    la = lb
    va(1:3) = vb(1:3)
    da = db
    go to 10
  end if
!
!  Try to find a starting edge based on adjacent high MDF values.
!  Select best of 4 possibilities based on max-min angle criterion.
!
  vave = val(1)
  m = 1

  do i = 2, n
    vave = vave + val(i)
    if ( val(m) < val(i) ) then
      m = i
    end if
  end do

  vave = vave / real ( n, kind = 8 )
  l = m
  i = m

30 continue

  i = i + 1
  if ( n < i ) then
    i = 1
  end if

  if ( vave < val(i) ) then
    m = i
    go to 30
  end if

  i = l

40 continue

  i = i - 1
  if ( i < 1 ) then
    i = n
  end if

  if ( vave < val(i) ) then
    l = i
    go to 40
  end if

  ll(1) = l
  ll(2) = l - 1
  if ( ll(2) < 1) ll(2) = n
  mm(1) = m
  mm(2) = m + 1

  if ( n < mm(2) ) then
    mm(2) = 1
  end if

  cosmin = 2.0D+00

  do i = 1, 2

    ii(1) = ll(i)
    s(1) = indv(ii(1))

    if ( 0 < s(1) ) then
      jj(1) = fvl(loc,s(1))
    else
      jj(1) = -s(1)
    end if

    do j = 1,2

      ii(2) = mm(j)

      if ( abs(ii(1) - ii(2)) <= 1 ) then
        cycle
      end if

      s(2) = indv(ii(2))

      if ( 0 < s(2) ) then
        jj(2) = fvl(loc,s(2))
      else
        jj(2) = -s(2)
      end if

      va(1:3) = vcl(1:3,jj(2)) - vcl(1:3,jj(1))

      da = va(1)**2 + va(2)**2 + va(3)**2
      cosmax = -1.0D+00

      do k = 1, 2

        if ( k == 2 ) then
          va(1:3) = -va(1:3)
        end if

        lb = ii(k) - 1
        if ( lb < 1 ) then
          lb = n
        end if
        lb = indv(lb)

        if ( 0 < lb ) then
          lb = fvl(loc,lb)
        else
          lb = -lb
        end if

        vb(1:3) = vcl(1:3,lb) - vcl(1:3,jj(k))

        db = vb(1)**2 + vb(2)**2 + vb(3)**2
        dotp = (va(1)*vb(1)+va(2)*vb(2)+va(3)*vb(3))/sqrt(da*db)
        cosmax = max ( cosmax, dotp )

        if ( s(k) < 0 ) then

          cosmax = max ( cosmax, -dotp )

        else

          lb = ii(k) + 1

          if ( n < lb ) then
            lb = 1
          end if

          lb = indv(lb)

          if ( 0 < lb ) then
            lb = fvl(loc,lb)
          else
            lb = -lb
          end if

          vb(1:3) = vcl(1:3,lb) - vcl(1:3,jj(k))

          db = vb(1)**2 + vb(2)**2 + vb(3)**2
          dotp = (va(1)*vb(1) + va(2)*vb(2) + va(3)*vb(3))/ sqrt(da*db)
          cosmax = max ( cosmax, dotp )

        end if

      end do

      if ( cosmax < cosmin ) then
        cosmin = cosmax
        l = ii(1)
        m = ii(2)
      end if

    end do

  end do

  if ( max ( mxcos, cos ( angacc ) ) < cosmin ) then
    return
  end if
!
!  For starting edge, try 3 different cut planes.
!
  neg = ( abs(facep(2,f)) /= p )

  if ( neg ) then
    sf = facep(3,f)
  else
    sf = facep(2,f)
  end if

  kvc = nvc
  j = indv(l)
  k = indv(m)

  if ( k < j .and. j < 0 ) then
    call i4_swap ( l, m )
    call i4_swap ( j, k )
  end if

  if ( 0 < k ) then

    lb = fvl(loc,k)

  else

    kvc = kvc + 1
    lb = -k

    vcl(1:3,kvc) = vcl(1:3,lb)

    lb = kvc
    m = m - 1
    if ( m < 1 ) then
      m = n
    end if

    k = indv(m)
    la = fvl(loc,k)
    m = fvl(loc,fvl(succ,k))

    if ( 0 < (m - la)*sf ) then
      k = fvl(edgc,k)
    else
      k = fvl(edga,k)
    end if

    cedge(2,0) = k

  end if

  if ( 0 < j ) then
    la = fvl(loc,j)
    if ( 0 < sf ) then
      j = fvl(pred,j)
    end if
  else
    kvc = kvc + 1
    la = -j

    vcl(1:3,kvc) = vcl(1:3,la)

    la = kvc
    l = l - 1
    if ( l < 1) l = n
    j = indv(l)
  end if

  cedge(1,0) = lb
  cedge(1,1) = la
  cedge(2,1) = -j

  edg(1:3) = vcl(1:3,lb) - vcl(1:3,la)

  da = sqrt(edg(1)**2 + edg(2)**2 + edg(3)**2)
  edg(1:3) = edg(1:3)/da

  enr(1) = edg(2)*nrml(3,f) - edg(3)*nrml(2,f)
  enr(2) = edg(3)*nrml(1,f) - edg(1)*nrml(3,f)
  enr(3) = edg(1)*nrml(2,f) - edg(2)*nrml(1,f)

  if ( neg ) then
    enr(1:3) = -enr(1:3)
  end if

  do i = 1, 3

    ang = pi * fract(i)

    if ( msglvl == 4 ) then
      write ( *,600) j,lb,la,f,p,ang*180.0D+00/ pi
    end if

    cdang(1) = ang
    cn = cos(ang)
    if ( neg) cn = -cn
    ce = sin(ang)

    nrmlc(1:3) = cn*nrml(1:3,f) + ce*enr(1:3)

    nrmlc(4) = nrmlc(1)*vcl(1,la) + nrmlc(2)*vcl(2,la) + &
      nrmlc(3)*vcl(3,la)

    call sepfac(p,cntr,nrmlc,angacc,mxcos,dtol,nvc,maxvc,vcl, &
      facep,nrml,fvl,eang,nce,cedge,cdang,aflag, ierr )

    if ( ierr /= 0 .or. aflag ) then
      return
    end if

  end do

  600 format (' sfc2mf: a,lb,la,f,p,ang =',5i5,f9.4)

  return
end
