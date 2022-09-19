subroutine prmdf3 ( ipolh, widp, nvc, vcl, nrml, ivrt, xivrt, ifac, xifac, &
  facval, edgval, vrtval, nfcev, nedev, nvrev, listev, infoev, htsiz, &
  maxedg, ht, edge, ierr )

!*****************************************************************************80
!
!! PRMDF3 does preprocessing for the mesh distribution function evaluation.
!
!  Discussion:
!
!    This routine is a preprocessing step for evaluating a mesh distribution
!    function in polyhedron IPOLH - the faces, edges, vertices for
!    which distances must be computed are determined.
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
!    Input, integer ( kind = 4 ) IPOLH, the index of polyhedron.
!
!    Input, real ( kind = 8 ) WIDP, the width of polyhedron IPOLH.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates in VCL array.
!
!    Input, integer ( kind = 4 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), the indices of face vertices in VCL,
!    ordered by face.
!
!    Input, integer ( kind = 4 ) XIVRT(1:*), the pointer to first vertex of each face
!    in IVRT; vertices of face K are IVRT(I) for I from XIVRT(K) to
!    XIVRT(K+1)-1.
!
!    Input, integer ( kind = 4 ) IFAC(1:*), the indices of polyhedron faces in FACEP,
!    ordered by polyhedron.
!
!    Input, integer ( kind = 4 ) XIFAC(1:*), the pointer to first face of each polyhedron
!    in IFAC; faces of polyhedron IPOLH are IFAC(I) for I from
!    XIFAC(IPOLH) to XIFAC(IPOLH+1)-1.
!
!    Input, real ( kind = 8 ) FACVAL(1:*), the value associated with each
!    face of decomposition.
!
!    Input, real ( kind = 8 ) EDGVAL(1:*), the value associated with each
!    edge of decomposition.
!
!    Input, real ( kind = 8 ) VRTVAL(1:*), the value associated with each
!    vertex of decomposition.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime
!    number which is greater than or equal to the number of vertices in
!    polyhedron.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array;
!    should be at least number of edges.
!
!    Output, integer ( kind = 4 ) NFCEV, the number of faces for which distances must
!    be evaluated.
!
!    Output, integer ( kind = 4 ) NEDEV, the number of edges for which distances must
!    be evaluated.
!
!    Output, integer ( kind = 4 ) NVREV, the number of vertices for which distances
!    must be evaluated.
!
!    Output, integer ( kind = 4 ) LISTEV(1:NFCEV+NEDEV+NVREV), the indices of above
!    faces, edges, vertices; first are NFCEV indices in FACVAL of faces,
!    then NEDEV indices in EDGVAL of edges, and NVREV indices
!    in VRTVAL of vertices.
!
!    Output, integer ( kind = 4 ) INFOEV(1:4,1:NFCEV+NEDEV), info for evaluation of distances
!    associated with faces and edges; first NFCEV entries are
!    plane equations with unit normal for above faces:
!    INFOEV(1,I)*X+INFOEV(2,I)*Y+INFOEV(3,I)*Z = INFOEV(4,I);
!    last NEDEV entries are displacement vector (DX,DY,DZ) =
!    INFOEV(1:3,J) of edge from vertex IVRT(LISTEV(J)) and
!    length of edge = INFOEV(4,J).
!
!    [Note: It is assumed there is enough space for above 2 arrays.]
!
!    Workspace, integer HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG), the hash table
!    and edge records used to determine entries of LISTEV.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg

  integer ( kind = 4 ) bptr
  integer ( kind = 4 ) e
  integer ( kind = 4 ) edge(4,maxedg)
  real    ( kind = 8 ) edgval(*)
  integer ( kind = 4 ) f
  real    ( kind = 8 ) facval(*)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac(*)
  integer ( kind = 4 ) ind
  real    ( kind = 8 ) infoev(4,*)
  integer ( kind = 4 ) ipolh
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) listev(*)
  integer ( kind = 4 ) nedev
  integer ( kind = 4 ) nfcev
  real    ( kind = 8 ) nrml(3,*)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvrev
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) vcl(3,*)
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) widp
  integer ( kind = 4 ) xifac(*)
  integer ( kind = 4 ) xivrt(*)

  ierr = 0
  k = 0

  do i = xifac(ipolh), xifac(ipolh+1)-1

    f = abs ( ifac(i) )

    if ( facval(f) < widp ) then

      k = k + 1
      listev(k) = f
      infoev(1,k) = nrml(1,f)
      infoev(2,k) = nrml(2,f)
      infoev(3,k) = nrml(3,f)
      j = ivrt(xivrt(f))
      infoev(4,k) = dot_product ( nrml(1:3,f), vcl(1:3,j) )

    end if

  end do

  nfcev = k

  hdfree = 0
  last = 0
  ht(0:htsiz-1) = 0

  do i = xifac(ipolh), xifac(ipolh+1)-1

    f = abs ( ifac(i) )
    l = xivrt(f+1) - 1
    e = l
    la = ivrt(l)

    do j = xivrt(f), l

      lb = ivrt(j)

      call edght ( la, lb, f, nvc, htsiz, maxedg, hdfree, last, ht, &
        edge, g, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      if ( 0 < g ) then
        if ( edgval(e) < min(facval(f),facval(g)) ) then
          k = k + 1
          listev(k) = e
          infoev(1:3,k) = vcl(1:3,lb) - vcl(1:3,la)
          infoev(4,k) = sqrt(infoev(1,k)**2 + infoev(2,k)**2 &
            + infoev(3,k)**2)
        end if
      end if

      la = lb
      e = j

    end do

  end do

  nedev = k - nfcev
!
!  HT, EDGE arrays are used below in a similar way to routine EDGHT
!  except that EDGE(1,*) is not used and no entries are deleted.
!
  last = 0
  ht(0:htsiz-1) = 0

  do i = xifac(ipolh), xifac(ipolh+1)-1

    f = abs ( ifac(i) )
    l = xivrt(f+1) - 1
    e = l

    do j = xivrt(f), l

      if ( edgval(e) < edgval(j) ) then
        g = e
      else
        g = j
      end if

      e = j
      lb = ivrt(j)
      ind = mod ( lb, htsiz )
      bptr = -1
      ptr = ht(ind)

60    continue

      if ( ptr /= 0 ) then

        if ( lb < edge(2,ptr) ) then
          go to 70
        else if ( edge(2,ptr) < lb ) then
          bptr = ptr
          ptr = edge(4,ptr)
          go to 60
        else
          if ( edgval(g) < edgval(edge(3,ptr)) ) then
            edge(3,ptr)=g
          end if
          cycle
        end if

      end if

70    continue

      last = last + 1

      if ( maxedg < last ) then
        ierr = 1
        return
      end if

      if ( bptr == -1 ) then
        ht(ind) = last
      else
        edge(4,bptr) = last
      end if

      edge(2,last) = lb
      edge(3,last) = g
      edge(4,last) = ptr

    end do

  end do

  do i = 1, last
    j = edge(2,i)
    e = edge(3,i)
    if ( vrtval(j) < edgval(e) ) then
      k = k + 1
      listev(k) = j
    end if
  end do

  nvrev = k - (nfcev + nedev)

  return
end
