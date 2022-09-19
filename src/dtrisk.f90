subroutine dtrisk ( k, npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, &
  bf_numac, nface, nsmplx, bf, fc, ht, iwk, wk, ierr )

!*****************************************************************************80
!
!! DTRISK constructs a Delaunay triangulation of vertices in KD.
!
!  Discussion:
!
!    This routine constructs a Delaunay triangulation of K-D vertices using
!    incremental approach and local transformations. Vertices are
!    first sorted in lexicographically increasing order, and
!    then are inserted one at a time from outside the convex hull.
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
!    Input, integer ( kind = 4 ) K, the dimension of points and triangulation.
!
!    Input, integer ( kind = 4 ) NPT, the number of K-D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT; a good choice is a
!    prime number which is about 1/8 * NFACE.
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) VM(1:NPT).  On input, indices of vertices of
!    VCL being triangulated.  On output, indices are permuted, so that
!    VCL(*,VM(1)), ... , VCL(*,VM(NPT)) are in lexicographic increasing order,
!    with possible slight reordering so first K+1 vertices
!    are not in same hyperplane
!
!    Output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array;
!    BF_NUM <= BF_MAX.
!
!    Output, integer ( kind = 4 ) NFC, the number of positions used in FC array;
!    NFC <= FC_MAX.
!
!    Output, integer ( kind = 4 ) BF_NUMAC, the number of boundary faces in triangulation;
!    BF_NUMAC <= BF_NUM.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation; NFACE <= NFC.
!
!    Output, integer ( kind = 4 ) NSMPLX, the number of simplices in triangulation.
!
!    Output, integer ( kind = 4 ) BF(1:K,1:BF_NUM), the array of boundary face records
!    containing pointers (indices) to FC; if FC(K+2,I) = -J < 0 and
!    FC(1:K,I) = ABC...G, then BF(1,J) points to other boundary face with
!    (K-2)-facet BC...G, BF(2,J) points to other boundary face
!    with facet AC...G, etc.; if BF(1,J) <= 0, record is not
!    used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:K+4,1:NFC), the array of face records which are in
!    linked lists in hash table with direct chaining. Fields are:
!    FC(1:K,*) - A,B,C,...,G with 1<=A<B<C<...<G<=NPT; indices
!    in VM of K vertices of face; if A<=0, record not used
!    (in linked list of avail records with indices <= NFC);
!    internal use: if B <= 0, face in queue, not in triang
!    FC(K+1:K+2,*) - D,E; indices in VM of (K+1)st vertex of 1
!    or 2 simplices with face ABC...G; if boundary face
!    then E < 0 and |E| is an index of BF array
!    FC(K+3,*) - HTLINK; pointer (index in FC) of next element
!    in linked list (or NULL = 0)
!    FC(K+4,*) - used internally for QLINK (link for queues or
!    stacks); pointer (index in FC) of next face in queue/
!    stack (or NULL = 0); QLINK = -1 indicates face is not
!    in any queue/stack, and is output value (for records
!    not in avail list), except:
!    FC(K+4,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC.
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and simplices of triangulation.
!
!    Workspace, integer IWK(1:6*K+2).
!
!    Workspace, real ( kind = 8 ) WK(1:K*K+2*K+1).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) alpha
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfp
  integer ( kind = 4 ) ctr
  integer ( kind = 4 ) fc(k+4,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) indf
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iwk(6*k+2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loc
  integer ( kind = 4 ) mat
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) mv
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) bf_numac
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) opsidk
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) s
  real    ( kind = 8 ) sum2
  integer ( kind = 4 ) top
  integer ( kind = 4 ) topnv
  real    ( kind = 8 ) vcl(k,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)
  real    ( kind = 8 ) wk(k*k+2*k+1)
  integer ( kind = 4 ) zpn
!
!  Permute elements of VM so that vertices are in lexicographic
!  order, and initialize data structures.
!
  ierr = 0

  km1 = k - 1
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  ctr = 1
  alpha = kp1
  mat = alpha + kp1
  ind = k
  indf = ind + kp1
  mv = indf + kp1
  loc = mv + kp1
  zpn = loc + k

  call dhpsrt ( k, npt, k, vcl, vm )

  call frsmpx ( k, .true., npt, vcl, vm, iwk, iwk(ind), wk(mat), ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( msglvl == 4 ) then
    write ( *,600) (vm(i),i=1,kp1)
    write ( *,610) (iwk(i),i=1,km1)
  end if

  do i = 1, k
    sum2 = vcl(i,vm(1))
    do j = 2, kp1
      sum2 = sum2 + vcl ( i,vm(j) )
    end do
    wk(i) = sum2 / real ( kp1, kind = 8 )
  end do

  ht(0:sizht-1) = 0

  hdavbf = 0
  hdavfc = 0
  bf_num = kp1
  nfc = kp1
  nsmplx = 1

  do i = 1, kp1
    j = 0
    do l = 1, kp1
      if ( l /= i ) then
        j = j + 1
        iwk(km1+j) = l
        bf(j,i) = l
      end if
    end do
    call htinsk ( k, i, iwk(ind), i, -i, npt, sizht, fc, ht )
  end do

  do j = 1, km1
    iwk(j) = iwk(j) + k - j
  end do

  s = 1
!
!  Insert I-th vertex into Delaunay triangle of first I-1 vertices.
!
  do i = kp2, npt

    vi = vm(i)

    if ( msglvl == 4 ) then
      write ( *,620) i,vi
    end if

    ip = i - 1
    if ( i == kp2 ) then
      ip = 2
    end if
    l = s

    do j = l, km1
      if ( i == iwk(j) ) then
        ip = j + 2
        s = j + 1
      end if
    end do
!
!  Form stacks of boundary faces involving vertex IP.
!  TOP is for stack of boundary faces to be tested for visibility.
!  FRONT is for stack of boundary faces visible from vertex I.
!  TOPNV is for stack of boundary faces not visible from I.
!
    if ( i == kp2 ) then

      top = 1
      do j = 1, k
        fc(kp4,j) = j + 1
      end do
      fc(kp4,kp1) = 0

    else

      if ( ip /= i - 1 ) then

        do l = 1, bf_num

          if ( bf(1,l) <= 0 ) then
            cycle
          end if

          do j = 1, k
            ptr = bf(j,l)
            do jj = 1, k
              if ( fc(jj,ptr) == ip ) then
                go to 120
              end if
            end do
          end do

        end do

120     continue

        bfi = ptr

      end if

        top = bfi
        back = bfi
        fc(kp4,back) = 0
        ptr = top

130     continue

        bfp = -fc(kp2,ptr)

        do j = 1, k
          if ( fc(j,ptr) /= ip ) then
            nbr = bf(j,bfp)
            if ( fc(kp4,nbr) == -1 ) then
              fc(kp4,back) = nbr
              back = nbr
              fc(kp4,back) = 0
            end if
          end if
        end do

        ptr = fc(kp4,ptr)

        if ( ptr /= 0 ) then
          go to 130
        end if

      end if
!
!  Find a boundary face visible from vertex I.
!
      front = 0
      topnv = 0

150   continue

      if ( top == 0 ) then
        go to 180
      end if

      ptr = top
      top = fc(kp4,ptr)

      do j = 1, k
        iwk(km1+j) = vm(fc(j,ptr))
      end do

      if ( opsidk(k,iwk(ind),vcl,.false.,wk(ctr),vcl(1,vi),wk(mat), &
        wk(alpha)) == 1 ) then

        front = ptr
        fc(kp4,ptr) = -1

170     continue

        if ( top == 0) then
          go to 180
        end if

        ptr = top
        top = fc(kp4,ptr)
        fc(kp4,ptr) = -1
        go to 170

      else

        fc(kp4,ptr) = topnv
        topnv = ptr

      end if

    go to 150

180 continue

    if ( front == 0 ) then
      ierr = 406
      return
    end if
!
!  Find remaining visible boundary faces, add new simplices with
!  vertex I, apply local transformation based on empty hypersphere crit.
!
    call vbfack ( k, vcl(1,vi), wk(ctr), vcl, vm, bf, fc, front, topnv, &
      iwk(ind), wk(mat), wk(alpha) )

    call nwsxou(k,i,npt,sizht,bf_num,nfc,bf_max,fc_max,bf,fc,ht,nsmplx, &
       hdavbf,hdavfc,front,back,bfi,iwk(ind), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call swaphs(k,i,npt,sizht,bf_num,nfc,bf_max,fc_max,vcl,vm,bf,fc,ht, &
      nsmplx,hdavbf,hdavfc,front,back,j,l,iwk(ind),iwk(indf), &
      iwk(mv),iwk(loc),iwk(zpn),wk(alpha),wk(mat), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( l /= 0 ) then
      bfi = l
    end if

  end do

  nface = nfc
  ptr = hdavfc

  do

    if ( ptr == 0 ) then
      exit
    end if

    nface = nface - 1
    ptr = -fc(1,ptr)

  end do

  bf_numac = bf_num
  ptr = hdavbf

  do

    if ( ptr == 0 ) then
      exit
    end if

    bf_numac = bf_numac - 1
    ptr = -bf(1,ptr)

  end do

  fc(kp4,1) = hdavbf
  fc(kp4,2) = hdavfc

  600 format (/1x,'dtrisk: first simplex: ',7i7)
  610 format (4x,'ishft(3:k+1)=',5i7)
  620 format (/1x,'step',i7,':   vertex i =',i7)

  return
end
