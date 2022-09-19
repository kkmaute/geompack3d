subroutine dsphdc ( nvc, nface, npolh, vcl, facep, nrml, fvl, eang, hfl, &
  pfl, htsiz, maxedg, edge, ht, ierr )

!*****************************************************************************80
!
!! DSPHDC initializes the polyhedral decomposition data structure.
!
!  Discussion:
!
!    This routine initializes the polyhedral decomposition data structure
!    where there are no holes on faces and no interior holes.  It is
!    assumed head vertex of each face is a strictly convex vertex.
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
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces in polyhedral decomposition.
!
!    Input, integer ( kind = 4 ) NPOLH, the number of polyhedra in decomposition.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1,1:NFACE+1), the head pointer to vertex indices
!    in FVL for each face; 1 = FACEP(1,1) < ... < FACEP(1,NFACE+1).
!
!    Input, integer ( kind = 4 ) FVL(1,1:*), the vertex indices; those for Ith face are in
!    FVL(1,J) for J = FACEP(1,I),...,FACEP(1,I+1)-1.
!
!    Input, integer ( kind = 4 ) HFL(1:NPOLH+1), the head pointer to face indices in PFL
!    for each polyhedron; 1 = HFL(1) < HFL(2) < ... < HFL(NPOLH+1).
!
!    Input, integer ( kind = 4 ) PFL(1,1:*), the signed face indices; those for Ith
!    polyhedron are in PFL(1,J) for J = HFL(I),...,HFL(I+1)-1; the face index
!    must be negated if the ordering of vertices for the face
!    in FVL is in CW order when viewed from outside Ith polyhedron.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime
!    number which is greater than or equal to NVC+2.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array;
!    should be at least maximum number of edges in a polyhedron of
!    decomposition.
!
!    Output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), FACEP(1,F) same as input;
!    FACEP(2,F) and FACEP(3,F) are signed indices of 2 polyhedra sharing face
!    F; if F is boundary face then FACEP(3,F) = 0; the sign of
!    the polyhedron index indicates whether face is oriented
!    counterclockwise (positive) or clockwise (negative) in
!    FVL when viewed from outside the polyhedron; if interior
!    face, 2 signs are different.
!
!    Output, real ( kind = 8 ) NRML(1:3,1:NFACE), the normals at faces;
!    NRML(*,F) is unit outward normal of face F with its vertices oriented
!    counterclockwise when viewed from outside polyhedron |FACEP(2,F)|.
!
!    Output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list where
!    NVERT = FACEP(1,NFACE+1)-1; 6 rows are for LOC, FACN, SUCC, PRED, EDGA,
!    EDGC; first 4 fields are the same as that used for the
!    convex polyhedron data structure (see routine DSCPH).
!    EDGA and EDGC give information about the edge UV where
!    V = FVL(SUCC,U). Let LU = FVL(LOC,U), LV = FVL(LOC,V),
!    and SF = +1 (-1) if face containing UV in polyhedron P is
!    oriented counterclockwise (CW) when viewed from outside P. Let WX be
!    edge corresponding to UV in the adjacent face of P, where
!    X = FVL(SUCC,W). If (LV-LU)*SF greater than 0, then FVL(EDGC,U) = W,
!    FVL(EDGA,W) = U, and EANG(U) is angle at UV between the
!    2 faces inside P; else FVL(EDGA,U) = W, FVL(EDGC,W) = U,
!    and EANG(W) is the edge angle. In other words, if P is
!    viewed from outside with edge UV directed upwards from
!    vertex with smaller LOC value to other vertex, then there
!    is a counterclockwise or CW rotation in P from face containing UV to
!    other face as indicated by EDGA or EDGC, respectively (A
!    for AntiCW, C for clockwise). If the counterclockwise
!    or clockwise rotation between 2 faces is exterior to the region,
!    then the EDGA or EDGC value is 0 and EANG value is -1.
!
!    Output, real ( kind = 8 ) EANG(1:NVERT), the angles at edges common
!    to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J) and is
!    determined by EDGC field.
!
!    Output, integer ( kind = 4 ) PFL(1:2,1:NPF); row 1 same as input and row 2 used
!    for link, where NPF = HFL(NPOLH+1)-1.
!
!    Workspace, integer HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG), the hash table
!    and edge records used to determine matching occurrences of polyhedron
!    edges by calling routine EDGHT.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npolh
  integer ( kind = 4 ) nvc

  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) ang
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) eang(*)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edge(4,maxedg)
  real    ( kind = 8 ) en(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,nface+1)
  integer ( kind = 4 ), parameter :: facn = 2
  logical              fflag
  integer ( kind = 4 ) fvl(6,*)
  integer ( kind = 4 ) g
  logical              gflag
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) hfl(npolh+1)
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lc
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nht
  real    ( kind = 8 ) nrml(3,nface)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,*)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pi2
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) sf
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,nvc)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  pi2 = 2.0D+00 * pi
  hdfree = 0
  last = 0
  nht = 0

  ht(0:htsiz-1) = 0

  do i = 1, nface
    facep(2,i) = 0
    facep(3,i) = 0
    k = facep(1,i)
    l = facep(1,i+1) - 1
    do j = k, l
      fvl(facn,j) = i
      fvl(succ,j) = j + 1
      fvl(pred,j) = j - 1
      fvl(edga,j) = 0
      fvl(edgc,j) = 0
      eang(j) = -1.0D+00
    end do
    fvl(succ,l) = k
    fvl(pred,k) = l
  end do

  do i = 1, npolh

    k = hfl(i)
    l = hfl(i+1) - 1

    do j = k, l
      pfl(2,j) = j + 1
      f = pfl(1,j)
      p = sign ( i, f )
      f = abs(f)
      if ( facep(2,f) == 0 ) then
        facep(2,f) = p
      else
        facep(3,f) = p
      end if
    end do

    pfl(2,l) = k

  end do

  do f = 1, nface
    if ( 0 < facep(2,f) * facep(3,f) ) then
      ierr = 321
      return
    end if
  end do
!
!  Compute normals for each face from orientation in FACEP(2,*).
!
  do f = 1, nface

    if ( 0 < facep(2,f) ) then
      ccw = succ
    else
      ccw = pred
    end if

    j = facep(1,f)
    lb = fvl(loc,j)
    lc = fvl(loc,fvl(ccw,j))
    la = fvl(loc,fvl(7-ccw,j))
    ab(1:3) = vcl(1:3,lb) - vcl(1:3,la)
    ac(1:3) = vcl(1:3,lc) - vcl(1:3,la)

    nrml(1,f) = ab(2) * ac(3) - ab(3) * ac(2)
    nrml(2,f) = ab(3) * ac(1) - ab(1) * ac(3)
    nrml(3,f) = ab(1) * ac(2) - ab(2) * ac(1)

    leng = sqrt ( sum ( nrml(1:3,f)**2 ) )

    if ( 0.0D+00 < leng ) then
      nrml(1:3,f) = nrml(1:3,f) / leng
    end if

  end do
!
!  Determine EDGA, EDGC fields and compute EANG values.
!
  do p = 1, npolh

    nht = 0

    do i = hfl(p), hfl(p+1)-1

      sf = pfl(1,i)
      f = abs(sf)

      do j = facep(1,f), facep(1,f+1)-1

        la = fvl(loc,j)
        lb = fvl(loc,fvl(succ,j))

        call edght ( la, lb, j, nvc, htsiz, maxedg, hdfree, last, ht, &
          edge, k, ierr )

        if ( ierr /= 0 ) then
          return
        end if

        if ( k <= 0 ) then

          nht = nht + 1

        else

          nht = nht - 1
          g = fvl(facn,k)
          dotp = dot_product ( nrml(1:3,f), nrml(1:3,g) )

          if ( 1.0D+00-tol < abs ( dotp ) ) then
            dotp = sign ( 1.0D+00, dotp )
          end if

          fflag = (abs(facep(2,f)) == p)
          gflag = (abs(facep(2,g)) == p)

          if ( fflag .neqv. gflag ) then
            dotp = -dotp
          end if

          ang = pi - acos ( dotp )
!
!  Determine whether edge angle is reflex.
!
          ab(1:3) = vcl(1:3,lb) - vcl(1:3,la)

          en(1) = nrml(2,f) * ab(3) - nrml(3,f) * ab(2)
          en(2) = nrml(3,f) * ab(1) - nrml(1,f) * ab(3)
          en(3) = nrml(1,f) * ab(2) - nrml(2,f) * ab(1)

          if ( fflag .neqv. ( 0 < sf ) ) then
            en(1:3) = -en(1:3)
          end if
!
!  AC = (midpoint of A and B) + EN - A
!
          do l = 1, 3
            ac(l) = 0.5D+00*(vcl(l,lb) - vcl(l,la)) + en(l)
          end do

          dotp = dot_product ( ac(1:3), nrml(1:3,g) )

          if ( .not. gflag ) then
            dotp = -dotp
          end if

          if ( 0.0D+00 < dotp ) then
            ang = pi2 - ang
          end if

          if ( 0 < (lb - la)*sf ) then
            fvl(edgc,j) = k
            fvl(edga,k) = j
            eang(j) = ang
          else
            fvl(edga,j) = k
            fvl(edgc,k) = j
            eang(k) = ang
          end if

        end if

      end do

    end do

    if ( nht /= 0 ) then
      ierr = 322
      return
    end if

  end do

  return
end
