subroutine dtrimk ( k, npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, &
  bf_numac, nface, nsmplx, bf, fc, ht, iwk, wk, ierr )

!*****************************************************************************80
!
!! DTRIMK constructs a Delaunay triangulation of points in KD.
!
!  Discussion:
!
!    This routine constructs the Delaunay triangulation of K-D vertices using
!    incremental approach and implicit local transformations, i.e.
!    simplices are first all deleted, then all added at each step.
!
!    Vertices are inserted one at a time in order given by VM array.
!    The initial simplices created due to a new vertex are obtained
!    by a walk through the triangulation until location of vertex is known.
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
!    Input/output, VM(1:NPT).  On input, the indices of vertices of VCL
!    being triangulated.  On output, the third to (K+1)st elements may
!    be swapped so that first K+1 vertices are not in same hyperplane.
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
!    containing pointers (indices) to FC; if FC(K+2,I) = -J < 0 and FC(1:K,I) =
!    ABC...G, then BF(1,J) points to other boundary face with
!    (K-2)-facet BC...G, BF(2,J) points to other boundary face
!    with facet AC...G, etc.; if BF(1,J) <= 0, record is not
!    used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:K+4,1:NFC), the array of face records which are
!    in linked lists in hash table with direct chaining. Fields are:
!    FC(1:K,*) - A,B,C,...,G with 1<=A<B<C<...<G<=NPT; indices
!      in VM of K vertices of face; if A<=0, record not used
!      (in linked list of avail records with indices <= NFC);
!      internal use: if B <= 0, face in queue, not in triang
!    FC(K+1:K+2,*) - D,E; indices in VM of (K+1)st vertex of 1
!      or 2 simplices with face ABC...G; if boundary face
!      then E < 0 and |E| is an index of BF array
!    FC(K+3,*) - HTLINK; pointer (index in FC) of next element
!      in linked list (or NULL = 0)
!    FC(K+4,*) - used internally for QLINK (link for queues or
!      stacks); pointer (index in FC) of next face in queue/
!      stack (or NULL = 0); QLINK = -1 indicates face is not
!      in any queue/stack, and is output value (for records
!      not in avail list), except:
!    FC(K+4,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and simplices of triangulation.
!
!    Workspace, integer IWK(3*K+1).
!
!    Workspace, real ( kind = 8 ) WK(K*K+2*K+1).
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
  logical              bflag
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
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(3*k+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loc
  integer ( kind = 4 ) mat
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) bf_numac
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) sum2
  integer ( kind = 4 ) top
  real    ( kind = 8 ) vcl(k,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)
  real    ( kind = 8 ) wk(k*k+2*k+1)
!
!  Find initial valid simplex, and initialize data structures.
!
  ierr = 0
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  ctr = 1
  alpha = kp1
  mat = alpha + kp1
  ind = 1
  indf = ind + k
  loc = indf + kp1

  call frsmpx ( k, .false., npt, vcl, vm, iwk, iwk(indf), wk(mat), ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( msglvl == 4 ) then
    write ( *,600) (vm(i),i=1,kp1)
    write ( *,610) (iwk(i),i=1,k-1)
  end if

  do i = 1, k
    sum2 = vcl(i,vm(1))
    do j = 2, kp1
      sum2 = sum2 + vcl(i,vm(j))
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
        iwk(j) = l
        bf(j,i) = l
      end if
    end do
    call htinsk(k,i,iwk(ind),i,-i,npt,sizht,fc,ht)
  end do

  ifac = kp1
!
!  Insert Ith vertex into Delaunay triangulation of first I-1 vertices.
!  Walk through triangulation to find location of vertex I, delete simplices
!  whose circumhypersphere contains I in interior, then add simplices
!  involving I by joining I to faces in stack.
!
  do i = kp2, npt

    vi = vm(i)

    if ( fc(kp2,ifac) == i-1 ) then
      ivrt = kp2
    else
      ivrt = kp1
    end if

    call walktk(k,vcl(1,vi),npt,sizht,nsmplx,vcl,vm,fc,ht,ifac, &
      ivrt,iwk(ind),iwk(indf),iwk(loc),wk(alpha),wk(mat), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( ivrt == 0 ) then

      if ( msglvl == 4 ) then
        write ( *,620) i,vi,'out'
      end if

      bflag = .true.
      top = 0
      front = ifac
      back = ifac
      call vbfack ( k, vcl(1,vi), wk(ctr), vcl, vm, bf, fc, front, 0, &
        iwk(ind), wk(mat), wk(alpha) )
      ptr = front

      do

        l = -fc(kp2,ptr)
        bf(1,l) = -hdavbf
        hdavbf = l
        fc(kp2,ptr) = i
        ptr = fc(kp4,ptr)

        if ( ptr == 0 ) then
          exit
        end if

      end do

    else if ( ivrt == 1 ) then

      if ( msglvl == 4 ) then
        write ( *,620) i,vi,'vert'
      end if

      ierr = 402

    else

      if ( msglvl == 4 ) then
        if ( kp1 <= ivrt ) then
          write ( *,620) i,vi,'in'
        else if ( ivrt == k ) then
          write ( *,620) i,vi,'face'
        else if ( ivrt < k ) then
          write ( *,620) i,vi,'edge'
        end if
      end if

      call lfcini(k,i,ifac,ivrt,iwk(indf),npt,sizht,bf,fc,ht, &
        nsmplx,hdavbf,hdavfc,bflag,front,back,top,iwk(ind), &
        iwk(loc), ierr )

    end if

    if ( ierr /= 0 ) then
      return
    end if

    call smpxda(k,i,npt,sizht,bf_num,nfc,bf_max,fc_max,vcl,vm,bf,fc,ht, &
      nsmplx,hdavbf,hdavfc,bflag,front,back,top,ifac,iwk(ind), &
      iwk(indf),wk(alpha),wk(mat), ierr )

    if ( ierr /= 0 ) then
      return
    end if

  end do

  nface = nfc
  ptr = hdavfc

  do while ( ptr /= 0 )
    nface = nface - 1
    ptr = -fc(1,ptr)
  end do

  bf_numac = bf_num
  ptr = hdavbf

  do while ( ptr /= 0 )
    bf_numac = bf_numac - 1
    ptr = -bf(1,ptr)
  end do

  fc(kp4,1) = hdavbf
  fc(kp4,2) = hdavfc

  600 format (/1x,'dtrimk: first simplex: ',7i7)
  610 format (4x,'iswap(3:k+1)=',5i7)
  620 format (/1x,'step',i7,':   vertex i =',i7,3x,a)

  return
end
