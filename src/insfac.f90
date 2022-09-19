subroutine insfac ( p, nrmlc, nce, cedge, cdang, nvc, nface, nvert, npolh, &
  npf, maxfp, maxfv, maxhf, maxpf, vcl, facep, factyp, nrml, fvl, eang, &
  hfl, pfl, ierr )

!*****************************************************************************80
!
!! INSFAC inserts a new cut face into a polyhedral decomposition.
!
!  Discussion:
!
!    This routine inserts a new face (cut face) in polyhedral decomposition
!    data structure. It is assumed that interior of face does not
!    intersect any other faces.
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
!    Input, real ( kind = 8 ) NRMLC(1:3), the unit normal vector of cut plane.
!
!    Input, integer ( kind = 4 ) NCE, the number of edges in cut face.
!
!    Input/output, integer ( kind = 4 ) CEDGE(1:2,0:NCE).  On input, CEDGE(1,I) is an
!    index of VCL, indices greater than NVC are new points;
!    CEDGE(2,I) = J indicates that edge of cut face ending at CEDGE(1,I) is
!    edge from J to FVL(SUCC,J) if J greater than 0; else if J < 0 then
!    edge of cut face ending at CEDGE(1,I) is a new edge and CEDGE(1,I) lies
!    on edge from -J to FVL(SUC,-J) and new edge lies in face FVL(FACN,-J);
!    CEDGE(2,I) always refers to an edge in the subpolyhedron
!    in negative half-space; CEDGE(1,NCE) = CEDGE(1,0);
!    CEDGE(2,0) is not input but is used temporarily.
!    On output, CEDGE(1:1,1:NCE) is updated to edges of cut face with respect to
!    positive half-space, CEDGE(2:2,1:NCE) has negative entries updated to
!    index of new edge.
!
!    Input, real ( kind = 8 ) CDANG(1:NCE), dihedral angles created by edges
!    of cut polygon in positive half-space; negative sign for angle I
!    indicates that face containing edge I is oriented CW in polyhedron P.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates
!    (excluding new ones)).
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
!    Input, integer ( kind = 4 ) MAXFP, the maximum size available for FACEP, FACTYP,
!    NRML arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXHF, the maximum size available for HFL array.
!
!    Input, integer ( kind = 4 ) MAXPF, the maximum size available for PFL array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NVC+?), the vertex coordinate list;
!    the new vertices to be inserted as indicated by CEDGE are after column NVC.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal
!    vectors for faces.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, integer ( kind = 4 ) EANG(1:NVERT), the edge angles.
!
!    Input/output, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices
!    in PFL for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices
!    for each polyhedron.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxhf
  integer ( kind = 4 ) nce

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ang
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  real    ( kind = 8 ) cdang(nce)
  integer ( kind = 4 ) cedge(2,0:nce)
  logical              docf
  logical              dof
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) fp
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) head
  integer ( kind = 4 ) hfl(maxhf)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  real    ( kind = 8 ) nrml(3,maxfp)
  real    ( kind = 8 ) nrmlc(3)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,maxpf)
  integer ( kind = 4 ) pind
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) sf
  integer ( kind = 4 ) sp
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) tailn
  integer ( kind = 4 ) tailp
  integer ( kind = 4 ) tailt
  real    ( kind = 8 ) vcl(3,*)
!
!  Insert new vertices and update CEDGE(2,*).
!
  ierr = 0

  do i = 0, nce-1

    if ( cedge(1,i) <= nvc ) then
      cycle
    end if

    if ( i == 0 ) then
      j = nce
    else
      j = i
    end if

    a = -cedge(2,j)
    call insvr3(a,nvc,nvert,maxfv,vcl,fvl,eang,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    if ( cdang(j) < 0.0D+00 ) then
      cedge(2,j) = -fvl(succ,a)
    end if

  end do
!
!  Insert new edges and update CEDGE(2,*).
!
  cedge(2,0) = cedge(2,nce)

  do i = 1, nce

    b = -cedge(2,i)

    if ( b < 0 ) then
      cycle
    end if

    f = fvl(facn,b)
    la = cedge(1,i-1)
    a = cedge(2,i-1)
!
!  This can only occur if I = 1.
!
    if ( a < 0 ) then

      a = -a

      if ( fvl(loc,a) == la ) then
        a = fvl(pred,a)
      else
        a = fvl(succ,a)
      end if

    end if

    do

      if ( fvl(loc,a) == la ) then
        a = fvl(pred,a)
        j = la - fvl(loc,a)
        sf = p
      else
        a = fvl(succ,a)
        j = fvl(loc,fvl(succ,a)) - la
        sf = -p
      end if

      if ( 0 < j * sf ) then
        a = fvl(edgc,a)
      else
        a = fvl(edga,a)
      end if

      fp = fvl(facn,a)

      if ( fp == f ) then
        exit
      end if

    end do

    if ( fvl(loc,a) == la ) then
      j = a
      a = fvl(succ,b)
      b = j
    else
      a = fvl(succ,a)
    end if

    call insed3(a,b,nface,nvert,npf,maxfp,maxfv,maxpf,facep, &
      factyp,nrml,fvl,eang,hfl,pfl,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    cedge(2,i) = a

  end do
!
!  Insert cut face into decomposition data structure. Subpolyhedron
!  in negative half space is numbered P, other is numbered NPOLH.
!
  nface = nface + 1
  npolh = npolh + 1
  npf = npf + 2

  if ( maxfv < nvert + nce ) then
    ierr = 15
    return
  else if ( maxfp < nface ) then
    ierr = 16
    return
  else if ( maxpf < npf ) then
    ierr = 17
    return
  else if ( maxhf < npolh ) then
    ierr = 18
    return
  end if

  nv = nvert
  facep(1,nface) = nvert + 1
  facep(2,nface) = p
  facep(3,nface) = -npolh
  factyp(nface) = 0
  nrml(1,nface) = nrmlc(1)
  nrml(2,nface) = nrmlc(2)
  nrml(3,nface) = nrmlc(3)

  do i = 0, nce-1
    nvert = nvert + 1
    fvl(loc,nvert) = cedge(1,i)
    fvl(facn,nvert) = nface
    fvl(succ,nvert) = nvert + 1
    fvl(pred,nvert) = nvert - 1
  end do

  fvl(succ,nvert) = facep(1,nface)
  fvl(pred,facep(1,nface)) = nvert
!
!  Set CEDGE(1,*) to edges of face in polyhedron NPOLH (after split). New
!  face is counterclockwise from outside new P which contains edges
!  of CEDGE(2,*).
!
  do i = 1, nce
    a = cedge(2,i)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))
    if ( 0.0D+00 < (lb - la) * cdang(i) ) then
      cedge(1,i) = fvl(edgc,a)
    else
      cedge(1,i) = fvl(edga,a)
    end if
  end do
!
!  Determine which faces of old P belong to new P or other new polyhedron,
!  and update HFL, PFL.  NP1 is used as a temporary polyhedron index.
!  Faces of old P are put into FACEP(2,*) field.
!  FACTYP(F) is set to -1 for double occurring faces.
!
  dof = .false.
  np1 = npolh + 1
  j = hfl(p)

60 continue

  f = abs ( pfl(1,j) )

  if ( abs ( facep(2,f) ) == abs ( facep(3,f) ) ) then
    facep(2,f) = np1
    facep(3,f) = -np1
    factyp(f) = -1
    dof = .true.
  else if ( abs ( facep(2,f) ) == p ) then
    facep(2,f) = sign ( np1, facep(2,f) )
  else
    i = facep(2,f)
    facep(2,f) = sign ( np1, facep(3,f) )
    facep(3,f) = i
    nrml(1,f) = -nrml(1,f)
    nrml(2,f) = -nrml(2,f)
    nrml(3,f) = -nrml(3,f)
  end if

  j = pfl(2,j)

  if ( j /= hfl(p)) go to 60

  do i = 1, nce

    j = 2
    f = fvl(facn,cedge(1,i))

    if ( factyp(f) == -1 ) then
      if ( fvl(loc,cedge(1,i)) /= fvl(loc,nv+i) ) then
        j = 3
      end if
    end if

    facep(j,f) = sign ( npolh, facep(j,f) )
    j = 2
    f = fvl(facn,cedge(2,i))

    if ( factyp(f) == -1 ) then
      if ( cdang(i) < 0.0D+00) j = 3
    end if

    facep(j,f) = sign ( p, facep(j,f) )

  end do

  pfl(1,npf-1) = nface
  pfl(1,npf) = -nface
  tailp = npf - 1
  tailn = npf
  tailt = hfl(p)
  ptr = pfl(2,tailt)
  pfl(2,tailt) = 0
  hfl(p) = npf - 1
  hfl(npolh) = npf

80 continue

  if ( ptr == 0) go to 110
  j = ptr
  sp = pfl(1,ptr)
  f = abs ( sp )
  ptr = pfl(2,ptr)

  if ( factyp(f) /= -1 .or. 0 < sp ) then
    k = 2
  else
    k = 3
  end if

  sf = facep(k,f)

  if ( abs ( sf ) == p ) then

    pfl(2,tailp) = j
    tailp = j

  else if ( abs ( sf ) == npolh ) then

    pfl(2,tailn) = j
    tailn = j

  else

    a = facep(1,f)
    la = fvl(loc,a)

90  continue

    b = fvl(succ,a)
    lb = fvl(loc,b)

    if ( 0 < (lb - la)*sp ) then
      c = fvl(edgc,a)
    else
      c = fvl(edga,a)
    end if

    g = fvl(facn,c)
    i = 2

    if ( factyp(g) == -1 ) then
      if ( fvl(loc,c) == la ) then
        if ( 0 < sp ) then
          i = 3
        end if
      else
        if ( sp < 0 ) then
          i = 3
        end if
      end if
    end if

    if ( abs ( facep(i,g) ) == p ) then
      pfl(2,tailp) = j
      tailp = j
      facep(k,f) = sign ( p, facep(k,f) )
      go to 100
    else if ( abs ( facep(i,g) ) == npolh ) then
      pfl(2,tailn) = j
      tailn = j
      facep(k,f) = sign ( npolh, facep(k,f) )
      go to 100
    end if

    a = b
    la = lb

    if ( a /= facep(1,f) ) then
      go to 90
    end if

    pfl(2,tailt) = j
    pfl(2,j) = 0
    tailt = j

100 continue

  end if

  go to 80

110 continue

  pfl(2,tailp) = npf - 1
  pfl(2,tailn) = npf
!
!  Check whether cut face occurs twice in same polyhedron.
!  Temporarily modify PRED field of edges of cut face.
!
  docf = .false.

  do i = 1, nce
    do j = 1, 2
      a = cedge(j,i)
      fvl(pred,a) = -fvl(pred,a)
    end do
  end do

  do i = 1, 2

    if ( i == 1 ) then
      head = npf - 1
      pind = p
    else
      head = npf
      pind = npolh
    end if

    ptr = pfl(2,head)

140 continue

    sf = pfl(1,ptr)
    f = abs ( sf )
    a = facep(1,f)
    la = fvl(loc,a)

150 continue

    b = fvl(succ,a)
    lb = fvl(loc,b)

    if ( 0 < fvl(pred,a) ) then

      if ( 0 < (lb - la)*sf ) then
        c = fvl(edgc,a)
      else
        c = fvl(edga,a)
      end if

      g = fvl(facn,c)
      k = 2

      if ( factyp(g) == -1 ) then
        if ( fvl(loc,c) == la ) then
          if ( 0 < sf ) then
            k = 3
          end if
        else
          if ( sf < 0) k = 3
        end if
      end if

      if ( abs ( facep(k,g) ) /= pind ) then
        docf = .true.
        go to 170
      end if

    end if

    a = b
    la = lb
    if ( a /= facep(1,f)) go to 150

    ptr = pfl(2,ptr)
    if ( ptr /= head ) go to 140

  end do
!
!  Reset PRED field of edges of cut face.
!
170 continue

  do i = 1, nce
    do j = 1, 2
      a = cedge(j,i)
      fvl(pred,a) = -fvl(pred,a)
    end do
  end do
!
!  Update EDGA, EDGC, and EANG fields.
!
  do i = 1, nce

    a = cedge(2,i)
    c = cedge(1,i)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))
    ang = abs ( cdang(i) )

    if ( 0.0D+00 < ( lb - la ) * cdang(i) ) then
      fvl(edgc,a) = nv + i
      fvl(edga,c) = nv + i
      fvl(edgc,nv+i) = c
      fvl(edga,nv+i) = a
      eang(a) = eang(a) - ang
      eang(nv+i) = ang
    else
      fvl(edgc,c) = nv + i
      fvl(edga,a) = nv + i
      fvl(edgc,nv+i) = a
      fvl(edga,nv+i) = c
      eang(nv+i) = eang(c) - ang
      eang(c) = ang
    end if

  end do
!
!  If DOF, reset FACTYP values of -1 to 0.
!
  if ( dof ) then

    do i = 1, 2

      if ( i == 1 ) then
        head = npf - 1
      else
        head = npf
      end if

      ptr = pfl(2,head)

      do

        f = abs ( pfl(1,ptr) )

        if ( factyp(f) == -1 ) then
          factyp(f) = 0
        end if

        ptr = pfl(2,ptr)

        if ( ptr == head ) then
          exit
        end if

      end do

    end do

  end if
!
!  If cut face is double occurring, set all faces to belong to
!  polyhedron P.
!
  if ( .not. docf) go to 240

  npolh = npolh - 1
  facep(3,nface) = -p
  tailp = pfl(2,npf-1)
  pfl(2,npf-1) = npf
  ptr = npf

230 continue

  sf = pfl(1,ptr)
  f = abs ( sf )

  if ( 0 < sf * facep(2,f) ) then
    facep(2,f) = sign ( p, sf )
  else
    facep(3,f) = -p
  end if

  if ( pfl(2,ptr) /= npf ) then
    ptr = pfl(2,ptr)
    go to 230
  else
    pfl(2,ptr) = tailp
  end if

240 continue

  if ( msglvl == 2 ) then
    write ( *,600) nce,facep(2,nface),facep(3,nface)
    do i = 1, nce
      la = fvl(loc,nv+i)
      write ( *,610) i,la,(vcl(j,la),j=1,3)
    end do
    write ( *, '(a)' ) ' '
  end if

  600 format (1x,'cut face: #edges, polyh(1:2) =',3i7)
  610 format (1x,2i7,3f15.7)

  return
end
