subroutine itris3 ( npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, nface, &
  ntetra, bf, fc, ht, ierr )

!*****************************************************************************80
!
!! ITRIS3 constructs an initial triangulation of 3D vertices.
!
!  Discussion:
!
!    This routine constructs an initial triangulation of 3D vertices by first
!    sorting them in lexicographically increasing (x,y,z) order and
!    then inserting 1 vertex at a time from outside the convex hull.
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
!    Input, integer ( kind = 4 ) NPT, number of 3D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, size of hash table HT; a good choice is a prime
!    number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.  On output, indices are permuted, so that VCL(*,VM(1)), ... ,
!    VCL(*,VM(NPT)) are in lexicographic increasing order,
!    with possible slight reordering so first 4 vertices are
!    non-coplanar.
!
!    Output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array;
!    BF_NUM <= BF_MAX.
!
!    Output, integer ( kind = 4 ) NFC, the number of positions used in FC array;
!    NFC <= FC_MAX.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation;
!    NFACE <= NFC.
!
!    Output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Output, integer ( kind = 4 ) BF(1:3,1:BF_NUM), the  array of boundary face records
!    containing pointers (indices) to FC; if FC(5,I) = -J < 0 and
!    FC(1:3,I) = ABC, then BF(1,J) points to other boundary face with edge BC,
!    BF(2,J) points to other boundary face with edge AC, and
!    BF(3,J) points to other boundary face with edge AB;
!    if BF(1,J) <= 0, record is not used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:NFC), the array of face records which are in
!    linked lists in hash table with direct chaining. Fields are:
!    FC(1:3,*) - integer A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
!    vertices of face; if A <= 0, record is not used (it is
!    in linked list of avail records with indices <= NFC);
!    internal use: if B <= 0, face in queue, not in triang
!    FC(4:5,*) - D,E; indices in VM of 4th vertex of 1 or 2
!    tetrahedra with face ABC; if ABC is boundary face
!    then E < 0 and |E| is an index of BF array
!    FC(6,*) - HTLINK; pointer (index in FC) of next element
!    in linked list (or NULL = 0)
!    FC(7,*) - used internally for QLINK (link for queues or
!    stacks); pointer (index in FC) of next face in queue/
!    stack (or NULL = 0); QLINK = -1 indicates face is not
!    in any queue/stack, and is output value (for records
!    not in avail list), except:
!    FC(7,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and tetrahedra of triangulation.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfi
  real    ( kind = 8 ) ctr(3)
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) op
  integer ( kind = 4 ) opside
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) top
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)
!
!  Permute elements of VM so that vertices are in lexicographic
!  order, and initialize data structures.
!
  call dhpsrt ( 3, npt, 3, vcl, vm )
!
!  Reorder points so that first four points are in general position.
!
  call frstet ( .true., npt, vcl, vm, i3, i4, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  do i = 1, 3
    ctr(i) = ( vcl(i,vm(1)) + vcl(i,vm(2)) + vcl(i,vm(3)) + &
      vcl(i,vm(4)) ) / 4.0D+00
  end do

  ht(0:sizht-1) = 0
  hdavbf = 0
  hdavfc = 0
  bf_num = 4
  nfc = 4
  ntetra = 1
  call htins ( 1, 1, 2, 3, 4, -1, npt, sizht, fc, ht )
  call htins ( 2, 1, 2, 4, 3, -2, npt, sizht, fc, ht )
  call htins ( 3, 1, 3, 4, 2, -3, npt, sizht, fc, ht )
  call htins ( 4, 2, 3, 4, 1, -4, npt, sizht, fc, ht )

  bf(1,1) = 4
  bf(2,1) = 3
  bf(3,1) = 2
  bf(1,2) = 4
  bf(2,2) = 3
  bf(3,2) = 1
  bf(1,3) = 4
  bf(2,3) = 2
  bf(3,3) = 1
  bf(1,4) = 3
  bf(2,4) = 2
  bf(3,4) = 1

  if ( msglvl == 4 ) then
    write ( *,600) (vm(i),i=1,4),i3,i4
  end if
!
!  Insert I-th vertex into triangulation of first I-1 vertices.
!
  do i = 5, npt

    vi = vm(i)
    if ( msglvl == 4 ) then
      write ( *,610) i,vi
    end if

    ip = i - 1

    if ( i == 5 ) then
      ip = 2
    end if

    if ( i == i3 + 2 ) then
      ip = 3
    end if

    if ( i == i4 + 1 ) then
      ip = 4
    end if
!
!  Form stacks of boundary faces involving vertex IP.
!  TOP is for stack of boundary faces to be tested for visibility.
!  FRONT is for stack of boundary faces visible from vertex I.
!  TOPNV is for stack of boundary faces not visible from I.
!
    front = 0
    topnv = 0

    if ( i == 5 ) then
      top = 4
      a = 3
      b = 2
      if ( ip == 2 ) then
        a = 2
      end if

      if ( ip <= 3 ) then
        b = 1
      end if

      fc(7,top) = a
      fc(7,a) = b
      fc(7,b) = 0

    else if ( ip == i - 1 ) then

      top = bfi
      fc(7,bfi) = 0
      b = fc(2,bfi)
      ptr = bf(1,-fc(5,bfi))

      do

        if ( fc(1,ptr) == b ) then
          b = fc(2,ptr)
          j = 1
        else
          b = fc(1,ptr)
          j = 2
        end if

        fc(7,ptr) = top
        top = ptr
        ptr = bf(j,-fc(5,ptr))

        if ( ptr == bfi ) then
          exit
        end if

      end do

    else

      do k = 1, bf_num

        if ( bf(1,k) <= 0 ) then
          cycle
        end if

        do e = 1, 3

          ptr = bf(e,k)

          if ( fc(1,ptr) == ip ) then
            b = fc(2,ptr)
            j = 3
            go to 60
          else if ( fc(2,ptr) == ip ) then
            b = fc(1,ptr)
            j = 3
            go to 60
          else if ( fc(3,ptr) == ip ) then
            b = fc(1,ptr)
            j = 2
            go to 60
          end if

        end do

      end do

60    continue

      bfi = ptr
      top = bfi
      fc(7,bfi) = 0
      ptr = bf(j,-fc(5,bfi))

70    continue

      if ( fc(1,ptr) == b ) then
        j = 1
        if ( fc(2,ptr) == ip ) then
          b = fc(3,ptr)
        else
          b = fc(2,ptr)
        end if
      else if ( fc(2,ptr) == b ) then
        j = 2
        if ( fc(1,ptr) == ip ) then
          b = fc(3,ptr)
        else
          b = fc(1,ptr)
        end if
      else
        j = 3
        if ( fc(1,ptr) == ip ) then
          b = fc(2,ptr)
        else
          b = fc(1,ptr)
        end if
      end if

      fc(7,ptr) = top
      top = ptr
      ptr = bf(j,-fc(5,ptr))
      if ( ptr /= bfi) go to 70

    end if
!
!  Find a boundary face visible from vertex I.
!
80  continue

    if ( top == 0) go to 110

    ptr = top
    top = fc(7,ptr)
    va = vm(fc(1,ptr))
    vb = vm(fc(2,ptr))
    vc = vm(fc(3,ptr))
    op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, vcl(1,vi) )

    if ( op == 2 ) then
      ierr = 301
      return
    end if

    if ( op == 1 ) then
      front = ptr
90    continue
      if ( top == 0) go to 110
      ptr = top
      top = fc(7,ptr)
      fc(7,ptr) = -1
      go to 90
    else
      fc(7,ptr) = topnv
      topnv = ptr
    end if

    go to 80

110 continue

    if ( front == 0 ) then
      ierr = 306
      return
    end if
!
!  Find remaining visible boundary faces, add new tetrahedron with vertex I.
!
    call vbfac ( vcl(1,vi), ctr, vcl, vm, bf, fc, front, topnv )

    if ( ierr /= 0 ) then
      return
    end if

    call nwthou ( i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, &
      ntetra, hdavbf, hdavfc, front, back, bfi, ierr )

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

  fc(7,1) = hdavbf
  fc(7,2) = hdavfc

  600 format (/1x,'itris3: first tetrahedron: ',4i7/4x,'i3, i4 =',2i7)
  610 format (/1x,'step',i7,':   vertex i =',i7)

  return
end
