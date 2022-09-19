subroutine dscph ( nvc, nface, vcl, hvl, nrml, fvl, eang, htsiz, maxedg, &
  edge, ht, ierr )

!*****************************************************************************80
!
!! DSCPH initalizes the convex polyhedron data structure.
!
!  Discussion:
!
!    This routine initializes a data structure for a convex polyhedron. It is
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
!    Input, integer ( kind = 4 ) NFACE, the number of faces in convex polyhedron.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NFACE+1), the head pointer to vertex indices in FVL
!    for each face; 1 = HVL(1) < HVL(2) < ... < HVL(NFACE+1).
!
!    Input, integer ( kind = 4 ) FVL(1,1:*), the vertex indices in counterclockwise
!    order when viewed from outside polyhedron; those for Ith face are in
!    FVL(1,J) for J = HVL(I),...,HVL(I+1)-1.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime
!    number which is greater than or equal to NVC+2.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array;
!    should be at least number of edges in polyhedron.
!
!    Output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit outward normals
!    of polyhedron faces.
!
!    Output, integer ( kind = 4 ) FVL(1:5,1:NVERT), EANG(1:NVERT), the face vertex list,
!    edge angles where NVERT = HVL(NFACE+1)-1, first row of FVL same as
!    input, and HVL(NFACE+1) not needed on output; contains
!    the 6 'arrays' LOC, FACN, SUCC, PRED, EDGV, EANG (first 5
!    are integer arrays, last is a double precision array);
!    the vertices of each face are stored in counterclockwise order (when
!    viewed from outside polyhedron) in doubly circular linked
!    list.  FVL(LOC,V) is location in VCL of the coordinates
!    of 'vertex' (index) V. EANG(V) is edge angle at edge
!    starting at vertex V. FVL(FACN,V) is face number (index
!    of HVL) of face containing V.  FVL(SUCC,V) [FVL(PRED,V)]
!    is index in FVL of the successor (predecessor) vertex of
!    vertex V.  FVL(EDGV,V) gives information about the edge
!    joining vertices V and its successor - it is equal to the
!    index in FVL of the successor vertex as represented in
!    the other face of the polyhedron sharing this edge, i.e.
!    FVL(EDGV,V) != FVL(SUCC,V), FVL(LOC,FVL(EDGV,V)) =
!    FVL(LOC,FVL(SUCC,V)), FVL(EDGV,FVL(EDGV,V)) = V.
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
  integer ( kind = 4 ) nvc

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) ang
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) eang(*)
  integer ( kind = 4 ) edge(4,maxedg)
  integer ( kind = 4 ), parameter :: edgv = 5
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(5,*)
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) hvl(nface+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) last
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nht
  real    ( kind = 8 ) nrml(3,nface)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,nvc)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  hdfree = 0
  last = 0
  nht = 0

  ht(0:htsiz-1) = 0

  do i = 1, nface

    k = hvl(i)
    l = hvl(i+1) - 1

    do j = k, l
      fvl(facn,j) = i
      fvl(succ,j) = j + 1
      fvl(pred,j) = j - 1
    end do

    fvl(succ,l) = k
    fvl(pred,k) = l

    do j = k, l

      call edght ( fvl(loc,j), fvl(loc,fvl(succ,j)), j, nvc, htsiz, &
        maxedg, hdfree, last, ht, edge, a, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      if ( 0 < a ) then
        fvl(edgv,j) = a
        fvl(edgv,a) = j
        nht = nht - 1
      else
        nht = nht + 1
      end if

    end do

  end do

  if ( nht /= 0 ) then
    ierr = 310
    return
  end if
!
!  Compute unit outward normals of faces.
!
  do f = 1, nface

    i = hvl(f)
    a = fvl(loc,fvl(pred,i))
    b = fvl(loc,i)
    c = fvl(loc,fvl(succ,i))
    ab(1:3) = vcl(1:3,b) - vcl(1:3,a)
    ac(1:3) = vcl(1:3,c) - vcl(1:3,a)

    nrml(1,f) = ab(2) * ac(3) - ab(3) * ac(2)
    nrml(2,f) = ab(3) * ac(1) - ab(1) * ac(3)
    nrml(3,f) = ab(1) * ac(2) - ab(2) * ac(1)

    leng = sqrt ( sum ( nrml(1:3,f)**2 ) )

    if ( 0.0D+00 < leng ) then
      nrml(1:3,f) = nrml(1:3,f)/leng
    end if

  end do
!
!  Compute angles at edges common to 2 faces.
!
  do f = 1, nface

    i = hvl(f)

70  continue

    a = fvl(edgv,i)
    j = fvl(facn,a)

    if ( f <= j ) then
      dotp = nrml(1,f)*nrml(1,j) + nrml(2,f)*nrml(2,j) + nrml(3,f)*nrml(3,j)
      if ( 1.0D+00 - tol < abs ( dotp ) ) then
        dotp = sign ( 1.0D+00, dotp )
      end if
      ang = pi - acos ( dotp )
      eang(i) = ang
      eang(a) = eang(i)
    end if

    i = fvl(succ,i)

    if ( i /= hvl(f) ) then
      go to 70
    end if

  end do

  return
end
