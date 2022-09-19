subroutine resedg ( u, angacc, rdacc, nvc, nface, nvert, npolh, npf, maxvc, &
  maxfp, maxfv, maxhf, maxpf, maxiw, maxwk, vcl, facep, factyp, nrml, fvl, &
  eang, hfl, pfl, rflag, iwk, wk, ierr )

!*****************************************************************************80
!
!! RESEDG resolves a reflex edge by a cut polygon.
!
!  Discussion:
!
!    This routine attempts to resolve a reflex edge by a cut polygon which
!    does not create small dihedral angles and does not pass near
!    any vertices not on cut plane.
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
!    Input, integer ( kind = 4 ) U, the index in FVL of reflex edge.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral
!    angle in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) RDACC, the minimum acceptable relative
!    distance between a cut plane and vertices not on plane.
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
!    Input/output, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions used
!    in HFL array.
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
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 5/3*NEDG where NEDG is number of edges in polyhedron
!    containing reflex edge.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array;
!    should be greater than or equal to NEDG.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex
!    coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit
!    normal vectors for faces.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Input/output, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices
!    in PFL for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices for
!    each polyhedron.
!
!    Output, logical RFLAG, TRUE iff reflex edge is resolved.
!
!    Workspace, integer IWK(1:3,1:MAXIW).
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
  integer ( kind = 4 ), parameter :: nangmx = 5

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ang(nangmx)
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) anghi
  real    ( kind = 8 ) anglo
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ca
  integer ( kind = 4 ) cb
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) ce
  real    ( kind = 8 ) cn
  real    ( kind = 8 ) da
  real    ( kind = 8 ) db
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dtol
  real    ( kind = 8 ) e(3)
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) enc(3)
  real    ( kind = 8 ) enl(3)
  real    ( kind = 8 ) enr(3)
  real    ( kind = 8 ) ep(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  logical              first
  integer ( kind = 4 ) fl
  integer ( kind = 4 ) fr
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) hfl(maxhf)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iwk(3,maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lc
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) luo
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lvo
  real    ( kind = 8 ) mindis
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nang
  integer ( kind = 4 ) nce
  integer ( kind = 4 ) nedg
  integer ( kind = 4 ) nedgc
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  real    ( kind = 8 ) nrml(3,maxfp)
  real    ( kind = 8 ) nrmlc(4)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,maxpf)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) rdacc
  logical              rflag
  integer ( kind = 4 ) sp
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) wk(maxwk)
!
!  Find faces FL, FR and polyhedron P containing reflex edge UV.
!  Do not resolve UV if FL or FR is a double-occurring face.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  rflag = .false.
  v = fvl(succ,u)
  lu = fvl(loc,u)
  lv = fvl(loc,v)
  fl = fvl(facn,u)
  fr = fvl(facn,fvl(edgc,u))
  if ( abs ( facep(2,fl) ) == abs ( facep(3,fl) ) ) then
    return
  end if

  if ( abs ( facep(2,fr) ) == abs ( facep(3,fr) ) ) then
    return
  end if

  if ( lu < lv ) then

    if ( 0 < facep(2,fl) ) then
      p = facep(2,fl)
    else
      p = facep(3,fl)
    end if

    luo = lu
    lvo = lv

  else

    if ( facep(2,fl) < 0 ) then
      p = -facep(2,fl)
    else
      p = -facep(3,fl)
    end if

    luo = lv
    lvo = lu

  end if

  if ( msglvl == 4 ) then
    write ( *,600) u,v,lu,lv,fl,fr,p
  end if
!
!  Determine edges in polyhedron P and average edge length.
!
  nedg = 0
  sum2 = 0.0D+00
  ptr = hfl(p)

10 continue

  f = pfl(1,ptr)

  if ( 0 < f ) then
    ccw = succ
  else
    ccw = pred
    f = -f
  end if

  a = facep(1,f)
  lb = fvl(loc,a)

20 continue

  la = lb
  b = fvl(ccw,a)
  lb = fvl(loc,b)

  if ( la < lb ) then

    leng = sqrt ( ( vcl(1,la) - vcl(1,lb) )**2 &
                + ( vcl(2,la) - vcl(2,lb) )**2 &
                + ( vcl(3,la) - vcl(3,lb) )**2 )

    sum2 = sum2 + leng

    if ( la /= luo .or. lb /= lvo ) then

      nedg = nedg + 1

      if ( maxiw < nedg ) then
        ierr = 6
        return
      end if

      if ( ccw == succ ) then
        iwk(1,nedg) = a
        iwk(2,nedg) = fvl(edgc,a)
      else
        iwk(1,nedg) = b
        iwk(2,nedg) = fvl(edgc,b)
      end if

    end if

  end if

  a = b
  if ( a /= facep(1,f) ) then
    go to 20
  end if

  ptr = pfl(2,ptr)

  if ( ptr /= hfl(p) ) then
    go to 10
  end if

  sum2 = sum2 / real ( nedg+1, kind = 8 )
  mindis = rdacc * sum2
  dtol = tol * sum2

  if ( 3 * maxiw < 5 * nedg ) then
    ierr = 6
    return
  else if ( maxwk < nedg ) then
    ierr = 7
    return
  end if
!
!  Compute unit vectors E, ENL, ENR which are directed edge UV and
!  edge normals to UV on faces FL and FR. E is unit normal vector
!  of plane containing ENL, ENR (it is normalization of ENL x ENR).
!
  e(1:3) = vcl(1:3,lvo) - vcl(1:3,luo)

  leng = sqrt(e(1)**2 + e(2)**2 + e(3)**2)

  e(1:3) = e(1:3) / leng

  enl(1) = nrml(2,fl)*e(3) - nrml(3,fl)*e(2)
  enl(2) = nrml(3,fl)*e(1) - nrml(1,fl)*e(3)
  enl(3) = nrml(1,fl)*e(2) - nrml(2,fl)*e(1)

  enr(1) = e(2) * nrml(3,fr) - e(3) * nrml(2,fr)
  enr(2) = e(3) * nrml(1,fr) - e(1) * nrml(3,fr)
  enr(3) = e(1) * nrml(2,fr) - e(2) * nrml(1,fr)

  if ( abs ( facep(2,fl) ) /= p ) then
    enl(1:3) = -enl(1:3)
  end if

  if ( abs ( facep(2,fr) ) /= p ) then
    enr(1:3) = -enr(1:3)
  end if

  k = 1
  if ( abs ( e(1) ) < abs ( e(2) ) ) then
    k = 2
  end if

  if ( abs ( e(k) ) < abs ( e(3) ) ) then
    k = 3
  end if
!
!  Find range of angles [ANGLO, ANGHI] that will resolve reflex edge,
!  and an ordered list of at most NANGMX angles determined by planes
!  containing reflex edge plus an adjacent edge (the bisection and
!  other angles may also be added to list). The list is obtained
!  by cycling through edges incident on U and V. counterclockwise = 1 (-1) if
!  going counterclockwise (CW) around LA (when viewed from
!  outside polyhedron P).
!
  ang(1) = eang(u) - pi
  anglo = max ( angacc, ang(1) )
  anghi = min ( pi, eang(u) - angacc)

  if ( ang(1) < anglo ) then
    nang = 0
  else
    nang = 2
    ang(2) = pi
  end if

  do i = 1, 2

    ccw = 1

    if ( i == 1 ) then
      la = lu
      b = fvl(pred,u)
      if ( lv < lu ) then
        ccw = -1
      end if
    else
      la = lv
      b = v
      if ( lu < lv ) then
        ccw = -1
      end if
    end if

    first = .true.
    f = fl

30  continue

    if ( .not. first ) then
      if ( fvl(loc,b) == la ) then
        b = fvl(pred,b)
      else
        b = fvl(succ,b)
      end if
    end if

    c = fvl(succ,b)

    if ( abs ( facep(2,f)) == abs ( facep(3,f) ) ) then

      if ( fvl(loc,b) == la ) then
        sp = -ccw*p
      else
        sp = ccw*p
       end if

    else if ( abs ( facep(2,f) ) == p ) then

      sp = facep(2,f)

    else

      sp = facep(3,f)

    end if

    if ( 0 < ( fvl(loc,c) - fvl(loc,b) ) * sp ) then
      b = fvl(edgc,b)
    else
      b = fvl(edga,b)
    end if

    f = fvl(facn,b)

    if ( f == fr ) then
      cycle
    end if

    if ( first ) then

      first = .false.

    else

      if ( fvl(loc,b) == la ) then
        c = fvl(succ,b)
      else
        c = b
      end if

      lc = fvl(loc,c)
      ep(1:3) = vcl(1:3,lc) - vcl(1:3,la)

      nrmlc(1) = ep(2) * e(3) - ep(3) * e(2)
      nrmlc(2) = ep(3) * e(1) - ep(1) * e(3)
      nrmlc(3) = ep(1) * e(2) - ep(2) * e(1)

      leng = sqrt ( sum ( nrmlc(1:3)**2 ) )

      if ( leng <= tol * &
        max ( abs ( ep(1) ), abs ( ep(2) ), abs ( ep(3) ) ) ) then
        go to 30
      end if

      nrmlc(1:3) = nrmlc(1:3) / leng

      enc(1) = e(2) * nrmlc(3) - e(3) * nrmlc(2)
      enc(2) = e(3) * nrmlc(1) - e(1) * nrmlc(3)
      enc(3) = e(1) * nrmlc(2) - e(2) * nrmlc(1)

      nrmlc(1) = enr(2)*enc(3) - enr(3)*enc(2)
      nrmlc(2) = enr(3)*enc(1) - enr(1)*enc(3)
      nrmlc(3) = enr(1)*enc(2) - enr(2)*enc(1)

      leng = sqrt ( sum ( nrmlc(1:3)**2 ) )

      if ( leng <= tol) then
        go to 30
      end if

      t = nrmlc(k)

      nrmlc(1) = enc(2)*enl(3) - enc(3)*enl(2)
      nrmlc(2) = enc(3)*enl(1) - enc(1)*enl(3)
      nrmlc(3) = enc(1)*enl(2) - enc(2)*enl(1)

      leng = sqrt ( sum ( nrmlc(1:3)**2 ) )

      if ( leng <= tol) then
        go to 30
      end if

      if ( nrmlc(k)*t < 0.0D+00) then
        go to 30
      end if

      dotp = enr(1)*enc(1) + enr(2)*enc(2) + enr(3)*enc(3)

      if ( e(k)*t < 0.0D+00 ) then
        dotp = -dotp
      end if

      if ( 1.0D+00 - tol < abs ( dotp ) ) then
        dotp = sign ( 1.0D+00, dotp )
      end if

      t = acos(dotp)

      if ( t < anglo .or. anghi < t ) then
        go to 30
      end if

      do j = 1, nang
        if ( abs ( ang(j) - t ) <= tol ) then
          go to 30
        end if
      end do

      nang = nang + 1
      ang(nang) = t
      if ( nangmx <= nang ) then
        go to 80
      end if

    end if

    go to 30

  end do

  t = eang(u) * 0.5D+00

  do j = 1, nang
    if ( abs ( ang(j) - t ) <= tol ) then
      go to 70
    end if
  end do

  nang = nang + 1
  ang(nang) = t

70 continue

  if ( nang <= nangmx - 2 .and. eang(u) <= 1.33D+00 * pi ) then

    nang = nang + 2
    ang(nang-1) = eang(u) * 0.25D+00
    ang(nang) = eang(u) * 0.75D+00

    if ( pi - ang(nang) <= 0.1D+00 ) then
      ang(nang-1) = eang(u) * 0.375D+00
      ang(nang) = eang(u) * 0.625D+00
    end if

    if ( ang(nang-1) < anglo .or. anghi < ang(nang-1) ) then
      ang(nang-1) = ang(nang)
      nang = nang - 1
    end if

    if ( ang(nang) < anglo .or. anghi < ang(nang) ) then
      nang = nang - 1
    end if

  else if ( nang < nangmx .and. eang(u) <= 1.4D+00 * pi ) then

    t = eang(u) * 0.375D+00

    if ( anglo <= t .and. t <= anghi ) then
      nang = nang + 1
      ang(nang) = t
    end if

  end if

80 continue

  if ( msglvl == 4 ) then
    write ( *,610) eang(u)*180.0D+00/ pi, &
      (ang(i)*180.0D+00 / pi,i=1,nang)
  end if
!
!  For each angle in ANG array, try to resolve reflex edge as
!  follows. Compute unit normal NRMLC and equation of cut plane; this
!  normal is outward with respect to subpolyhedron containing face FL; NRMLC(4)
!  is right hand side constant term of plane equation.
!  Then determine which edges of polyhedron intersect cut plane;
!  reject plane if there is a vertex within distance MINDIS of plane
!  which does not lie on plane.
!
  do i = 1, nang

    cn = cos(ang(i))
    if ( abs ( facep(2,fr) ) /= p ) then
      cn = -cn
    end if

    ce = sin(ang(i))

    nrmlc(1:3) = cn * nrml(1:3,fr) + ce*enr(1:3)

    nrmlc(4) = nrmlc(1)*vcl(1,lu) + nrmlc(2)*vcl(2,lu) + &
      nrmlc(3)*vcl(3,lu)

    nedgc = nedg
    j = 1

90  continue

    la = fvl(loc,iwk(1,j))
    da = nrmlc(1)*vcl(1,la) + nrmlc(2)*vcl(2,la) + &
      nrmlc(3)*vcl(3,la) - nrmlc(4)

    if ( abs(da) <= dtol ) then
      ca = 0
    else if ( abs(da) <= mindis ) then
      cycle
    else if ( da < 0.0D+00 ) then
      ca = -1
    else
      ca = 1
    end if

    b = fvl(succ,iwk(1,j))
    lb = fvl(loc,b)
    db = nrmlc(1)*vcl(1,lb) + nrmlc(2)*vcl(2,lb) + &
      nrmlc(3)*vcl(3,lb) - nrmlc(4)

    if ( abs(db) <= dtol ) then
      cb = 0
    else if ( abs(db) <= mindis ) then
      cycle
    else if ( db < 0.0D+00 ) then
      cb = -1
    else
      cb = 1
    end if

    if ( 0 < ca * cb ) then
      do k = 1, 2
        c = iwk(k,j)
        iwk(k,j) = iwk(k,nedgc)
        iwk(k,nedgc) = c
      end do
      nedgc = nedgc - 1
    else
      iwk(3,j) = (ca+2)*10 + (cb+2)
      j = j + 1
    end if

    if ( j <= nedgc ) then
      go to 90
    end if

    iwk(1,nedg+1) = lvo
    iwk(3,nedg+1) = luo
    iwk(1,nedg+2) = u
    wk(1) = ang(i)

    if ( msglvl == 4 ) then
      write ( *,620) i,nedgc
    end if

    call cutfac(p,nrmlc,angacc,dtol,nvc,maxvc,vcl,facep,nrml, &
      fvl,eang,nedgc,iwk,nce,iwk(1,nedg+1),wk,rflag,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    if ( rflag ) then
      call insfac(p,nrmlc,nce,iwk(1,nedg+1),wk,nvc,nface,nvert, &
        npolh,npf,maxfp,maxfv,maxhf,maxpf,vcl,facep,factyp,nrml, &
        fvl,eang,hfl,pfl,ierr)
      return
    end if

  end do

  600 format (' resedg: u,v,lu,lv,fl,fr,p =',7i5)
  610 format (4x,'candidate angles (dihedral angle = ',f12.5,')'/ &
     2x,f14.5,4f12.5)
  620 format (4x,'trying angle',i3,'   nedgc=',i5)

  return
end
