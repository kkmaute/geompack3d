subroutine smpxda ( k, i, npt, sizht, bf_num, nfc, bf_max, fc_max, vcl, vm, bf, &
  fc, ht, nsmplx, hdavbf, hdavfc, bflag, front, back, top, ifac, ind, indf, &
  center, mat, ierr )

!*****************************************************************************80
!
!! SMPXDA deletes simplices whose circumhypersphere contains a vertex.
!
!  Discussion:
!
!    This routine deletes simplices whose circumhypersphere contains a vertex
!    I in the interior, then adds simplices involving I by joining I to
!    faces in the stack, where I is index of new vertex added to triang.
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
!    Input, integer ( kind = 4 ) K, the dimension of triangulation.
!
!    Input, integer ( kind = 4 ) I, the local index of next vertex inserted in triangulation;
!    it is assumed I is largest index so far.
!
!    Input, integer ( kind = 4 ) NPT, the number of K-D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) BF(1:K,1:BF_MAX), the array of boundary face records;
!    see DTRISK.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:FC_MAX), the array of face records; see
!    routine DTRISK.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NSMPLX, the number of simplices in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVBF, the head pointer to available BF records.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input, logical BFLAG, TRUE iff vertex I is on boundary of triangulation.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, indices of front and back of queue
!    of interior faces that may form a new simplex with vertex I.
!
!    Input/output, integer ( kind = 4 ) TOP, index of top of stack of faces that form a
!    new simplex with vertex I.
!
!    Output, integer ( kind = 4 ) IFAC, the index of face containing vertex I.
!
!    Workspace, integer IND(1:K), the local vertex indices or pivot indices.
!
!    Workspace, integer INDF(1:K+1), the indices in VCL.
!
!    Workspace, real ( kind = 8 ) CENTER(1:K), the hypersphere center or
!    hyperplane normal.
!
!    Workspace, real ( kind = 8 ) MAT(1:K,1:K), the matrix used for solving
!    system of linear equations.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,bf_max)
  logical              bflag
  integer ( kind = 4 ) bfp
  real    ( kind = 8 ) center(k)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) fc(k+4,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) indf(k+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) l
  real    ( kind = 8 ) mat(k,k)
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) opsidk
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) rsq
  integer ( kind = 4 ) top
  integer ( kind = 4 ) topn
  real    ( kind = 8 ) vcl(k,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  vi = vm(i)
  topn = 0

10 continue

  if ( front == 0) go to 60

  pos = front
  front = fc(kp4,pos)

  if ( fc(2,pos) == 0 ) then
    fc(1,pos) = -hdavfc
    hdavfc = pos
    go to 10
  end if

  indf(1:k) = vm(fc(1:k,pos))

  if ( fc(kp1,pos) == i ) then
    d = fc(kp2,pos)
  else
    d = fc(kp1,pos)
  end if

  indf(kp1) = vm(d)
  call ccsphk ( k, .true., indf, vcl, vcl(1,vi), center, rsq, in, mat, ind )

  if ( in < 1 ) then
    fc(kp4,pos) = top
    top = pos
    go to 10
  end if
!
!  Delete simplex and process other faces of simplex.
!
  nsmplx = nsmplx - 1

  if ( msglvl == 4 ) then
    write ( *,600) (fc(ii,pos),ii=1,k),d
  end if

  do j = 1,k

    ind(1:k) = fc(1:k,pos)
    a = ind(j)
    ind(j) = d
    ptr = htsrck(k,ind,npt,sizht,fc,ht)
    if ( ptr <= 0) go to 160

    if ( 0 <= fc(kp4,ptr) ) then

      call htdelk ( k, ptr, npt, sizht, fc, ht )
      fc(2,ptr) = 0

    else if ( fc(kp2,ptr) < 0 ) then

      if ( bflag ) then

        indf(1:k) = vm(fc(1:k,ptr))

        if ( opsidk(k,indf,vcl,.true.,vcl(1,vi),vcl(1,vi),mat,center) == 0 ) then
          bfp = -fc(kp2,ptr)
          bf(1,bfp) = -hdavbf
          hdavbf = bfp
          call htdelk(k,ptr,npt,sizht,fc,ht)
          fc(1,ptr) = -hdavfc
          hdavfc = ptr
          cycle
        end if

      end if

      fc(kp1,ptr) = i
      fc(kp4,ptr) = top
      top = ptr

    else

      if ( fc(kp1,ptr) == a ) then
        fc(kp1,ptr) = i
      else
        fc(kp2,ptr) = i
      end if

      if ( front == 0 ) then
        front = ptr
      else
        fc(kp4,back) = ptr
      end if

      back = ptr
      fc(kp4,back) = 0

    end if

  end do

  call htdelk(k,pos,npt,sizht,fc,ht)
  fc(1,pos) = -hdavfc
  hdavfc = pos
  go to 10
!
!  For faces in stack TOP, form new simplices with vertex I.
!  Then set BF fields for new boundary faces if BFLAG = TRUE.
!
60 continue

  pos = top
  top = fc(kp4,pos)
  fc(kp4,pos) = -1
  nsmplx = nsmplx + 1

  if ( msglvl == 4 ) then
    write ( *,610) (fc(ii,pos),ii=1,k),i
  end if

  do j = 1, k

    ind(1:k) = fc(1:k,pos)
    a = ind(j)
    ind(j) = i
    ptr = htsrck(k,ind,npt,sizht,fc,ht)

    if ( 0 < ptr ) then

      fc(kp2,ptr) = a

    else

      call availk(k,hdavfc,nfc,fc_max,fc,ptr,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      call htinsk(k,ptr,ind,a,0,npt,sizht,fc,ht)

      if ( bflag ) then
        fc(kp4,ptr) = topn
        topn = ptr
      end if

    end if

  end do

  if ( top /= 0) go to 60

  ifac = ptr

  if ( .not. bflag ) then
    return
  end if

90 continue

  if ( topn == 0) go to 110

  pos = topn
  topn = fc(kp4,pos)

  if ( fc(kp2,pos) /= 0 ) then
    fc(kp4,pos) = -1
    go to 90
  end if

  if ( hdavbf /= 0 ) then

    bfp = hdavbf
    hdavbf = -bf(1,hdavbf)

  else

    if ( bf_max <= bf_num ) then
      ierr = 23
      return
    else
      bf_num = bf_num + 1
      bfp = bf_num
    end if

  end if

  fc(kp2,pos) = -bfp
  fc(kp4,pos) = top
  top = pos
  bf(1:k,bfp) = 0
  go to 90

110 continue

  pos = top
  top = fc(kp4,pos)
  fc(kp4,pos) = -1
  d = fc(kp1,pos)
  bfp = -fc(kp2,pos)

  do j = 1,k

    if ( bf(j,bfp) /= 0 ) then
      cycle
    end if

    a = fc(j,pos)
    b = d

120 continue

    ind(1:k) = fc(1:k,pos)
    ind(j) = b
    ptr = htsrck(k,ind,npt,sizht,fc,ht)

    if ( ptr <= 0) go to 160

    l = fc(kp2,ptr)

    if ( 0 < l ) then

      if ( l == a ) then
        a = b
        b = fc(kp1,ptr)
      else
        a = b
        b = fc(kp2,ptr)
      end if

      go to 120

    end if

    bf(j,bfp) = ptr
    ii = 1

140 continue

    if ( fc(ii,ptr) /= b ) then
      ii = ii + 1
      go to 140
    end if

    bf(ii,-l) = pos

  end do

  if ( top /= 0) go to 110
  return

160 continue

  ierr = 400

  600 format (1x,'deleted simplex:',9i7)
  610 format (1x,'added simplex:',9i7)

  return
end
