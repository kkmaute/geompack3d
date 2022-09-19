subroutine dsconv ( p, headp, facep, nrml, fvl, eang, pfl, ncface, ncvert, &
  chvl, cnrml, cfvl, ceang )

!*****************************************************************************80
!
!! DSCONV converts the representation of a convex polyhedron.
!
!  Discussion:
!
!    This routine converts the representation of a convex polyhedron in
!    the polyhedral decomposition data structure to the data
!    structure for a single convex polyhedron.
!
!    It is assumed upper bounds for NCVERT and NCFACE are
!    computed before calling this routine, and there is enough
!    space in CHVL, CNRML, CFVL, CEANG.
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
!    Input, integer ( kind = 4 ) P, the polyhedron index.
!
!    Input, integer ( kind = 4 ) HEADP, the head pointer to face indices in PFL for
!    polyhedron P.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list: row 1 is head
!    pointer, rows 2 and 3 are signed polyhedron indices.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces;
!    outward normal corresponds to counterclockwise traversal
!    of face from polyhedron with index |FACEP(2,F)|.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list; see routine DSPHDC,
!
!    Input, real ( kind = 8 ) EANG(1:*), angles at edges common to 2 faces
!    in a polyhedron; EANG(J) corresponds to FVL(*,J), determined by EDGC field
!    corresponds to FVL(*,J) and is determined by EDGC field.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:*), the list of signed face indices for each
!    polyhedron; row 2 used for link.
!
!    Output, integer ( kind = 4 ) NCFACE, the number of faces in convex polyhedron.
!
!    Output, integer ( kind = 4 ) NCVERT, the size of CFVL, CEANG arrays; 2 * number
!    of edges of polyhedron.
!
!    Output, integer ( kind = 4 ) CHVL(1:NCFACE), the head vertex list.
!
!    Output, real ( kind = 8 ) CNRML(1:3,1:NCFACE), the unit outward
!    normals of faces.
!
!    Output, integer ( kind = 4 ) CFVL(1:5,1:NCVERT), face vertex list; see routine DSCPH.
!
!    Output, real ( kind = 8 ) CEANG(1:NCVERT), the angles at edges common
!    to 2 faces; CEANG(I) corresponds to CFVL(*,I).
!
  implicit none

  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) ceang(*)
  integer ( kind = 4 ) cfvl(5,*)
  integer ( kind = 4 ) chvl(*)
  real    ( kind = 8 ) cnrml(3,*)
  real    ( kind = 8 ) eang(*)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ), parameter :: edgv = 5
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(6,*)
  integer ( kind = 4 ) headp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) kv
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lk
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) ncface
  integer ( kind = 4 ) ncvert
  real    ( kind = 8 ) nrml(3,*)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,*)
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) sf
  integer ( kind = 4 ), parameter :: succ = 3

  kf = 0
  kv = 0
  i = headp

10 continue

  sf = pfl(1,i)

  if ( 0 < sf ) then
    ccw = succ
  else
    ccw = pred
  end if

  f = abs(sf)
  kf = kf + 1
  chvl(kf) = kv + 1
  j = facep(1,f)

20 continue

  kv = kv + 1
  lj = fvl(loc,j)
  k = fvl(ccw,j)
  lk = fvl(loc,k)
  cfvl(loc,kv) = lj
  cfvl(facn,kv) = kf
  cfvl(succ,kv) = kv + 1
  cfvl(pred,kv) = kv - 1
  cfvl(edgv,kv) = 0

  if ( ccw == succ ) then
    fvl(facn,j) = kv
    r = j
  else
    fvl(facn,k) = kv
    s = lj
    lj = lk
    lk = s
    r = k
  end if

  if ( 0 < (lk - lj)*sf ) then
    ceang(kv) = eang(r)
  else
    ceang(kv) = eang(fvl(edga,r))
  end if

  j = fvl(ccw,j)

  if ( j /= facep(1,f) ) then
    go to 20
  end if

  cfvl(succ,kv) = chvl(kf)
  cfvl(pred,chvl(kf)) = kv

  if ( abs(facep(2,f)) == p ) then
    cnrml(1,kf) = nrml(1,f)
    cnrml(2,kf) = nrml(2,f)
    cnrml(3,kf) = nrml(3,f)
  else
    cnrml(1,kf) = -nrml(1,f)
    cnrml(2,kf) = -nrml(2,f)
    cnrml(3,kf) = -nrml(3,f)
  end if

  i = pfl(2,i)

  if ( i /= headp ) then
    go to 10
  end if

  ncface = kf
  ncvert = kv
!
!  Set CFVL(EDGV,*) field and reset FVL(FACN,*) field.
!
  i = headp

30 continue

  sf = pfl(1,i)
  f = abs(sf)
  j = facep(1,f)

40 continue

  r = fvl(facn,j)
  fvl(facn,j) = f

  if ( cfvl(edgv,r) == 0 ) then

    lj = fvl(loc,j)
    lk = fvl(loc,fvl(succ,j))

    if ( 0 < (lk - lj)*sf ) then
      s = fvl(facn,fvl(edgc,j))
    else
      s = fvl(facn,fvl(edga,j))
    end if

    cfvl(edgv,r) = s
    cfvl(edgv,s) = r
  end if

  j = fvl(succ,j)

  if ( j /= facep(1,f) ) then
    go to 40
  end if

  i = pfl(2,i)

  if ( i /= headp ) then
    go to 30
  end if

  return
end
