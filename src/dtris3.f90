subroutine dtris3 ( npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, nface, &
  ntetra, bf, fc, ht, ierr )

!*****************************************************************************80
!
!! DTRIS3 constructs a Delaunay triangulation of vertices in 3D.
!
!  Discussion:
!
!    This routine constructs a Delaunay triangulation of 3D vertices using
!    an incremental approach and local transformations.  Vertices are
!    first sorted in lexicographically increasing (x,y,z) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Modified:
!
!    02 September 2005
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
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of the hash table HT; a good choice is
!    a prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!    This needs to be at least as big as the number of boundary faces.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!    This needs to be at least as big as the number of faces.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NPT), the vertex coordinates.
!    In the general case, VCL may contain the coordinates for more
!    than NPT vertices, and the VM array is used to select them.
!
!    Input/output, integer ( kind = 4 ) VM(1:NPT), the vertices of VCL to be triangulated.
!    On output, these indices are permuted, so that VCL(*,VM(1)), ... ,
!    VCL(*,VM(NPT)) are in lexicographic increasing order,
!    with possible slight reordering so first 4 vertices are
!    non-coplanar.  Typically, the input value of VM might be 1 through
!    NPT.
!
!    Output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array;
!    BF_NUM <= BF_MAX.
!
!    Output, integer ( kind = 4 ) NFC, the number of positions used in FC array;
!    NFC <= FC_MAX.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation; NFACE <= NFC.
!
!    Output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in the triangulation.
!
!    Output, integer ( kind = 4 ) BF(1:3,1:BF_NUM), boundary face records containing pointers
!    (indices) to FC; if FC(5,I) = -J < 0 and FC(1:3,I) = ABC,
!    then BF(1,J) points to other boundary face with edge BC,
!    BF(2,J) points to other boundary face with edge AC, and
!    BF(3,J) points to other boundary face with edge AB;
!    if BF(1,J) <= 0, record is not used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:NFC), face records which are in linked lists
!    in hash table with direct chaining. Fields are:
!    FC(1:3,*) - A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
!    vertices of face; if A <= 0, record is not used (it is
!    in linked list of available records with indices <= NFC);
!    internal use: if B <= 0, face in queue, not in triangulation.
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
!    FC(7,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC.
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), a hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and tetrahedra of the triangulation.
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

  ierr = 0
!
!  Permute elements of VM so that vertices are in lexicographic order.
!
  call dhpsrt ( 3, npt, 3, vcl, vm )
!
!  Reorder points so that first four points are in general position.
!
  call frstet ( .true., npt, vcl, vm, i3, i4, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTRIS3 - Error!'
    write ( *, '(a,i6)' ) '  FRSTET returned IERR = ', ierr
    return
  end if
!
!  Initialize data structures.
!
  do i = 1, 3
    ctr(i) = sum ( vcl(i,vm(1:4)) ) / 4.0D+00
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

  bf(1:3,1) = (/ 4, 3, 2 /)
  bf(1:3,2) = (/ 4, 3, 1 /)
  bf(1:3,3) = (/ 4, 2, 1 /)
  bf(1:3,4) = (/ 3, 2, 1 /)

  if ( msglvl == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTRIS3:'
    write ( *, '(a)' ) '  First tetrahedron:'
    write ( *, '(2x,i6,2x,i6,2x,i6,2x,i6)' ) vm(1:4)
    write ( *, '(a,i6,a,i6)' ) '  I3 = ', i3, '  I4 = ', i4
  end if
!
!  Insert the I-th vertex into Delaunay triangle of first I-1 vertices.
!
  do i = 5, npt

    vi = vm(i)

    if ( msglvl == 4 ) then
      write ( *, '(a,i6,a,i6)' ) '  Step: ', i, '  Vertex: ', vi
    end if

    if ( i == 5 ) then
      ip = 2
    else
      ip = i - 1
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

      if ( ip == 2 ) then
        a = 2
      else
        a = 3
      end if

      if ( ip <= 3 ) then
        b = 1
      else
        b = 2
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

      j = 0

      do k = 1, bf_num

        if ( bf(1,k) <= 0 ) then
          cycle
        end if

        do e = 1, 3

          ptr = bf(e,k)

          if ( fc(1,ptr) == ip ) then
            b = fc(2,ptr)
            j = 3
            exit
          else if ( fc(2,ptr) == ip ) then
            b = fc(1,ptr)
            j = 3
            exit
          else if ( fc(3,ptr) == ip ) then
            b = fc(1,ptr)
            j = 2
            exit
          end if

        end do

        if ( j /= 0 ) then
          exit
        end if

      end do

      bfi = ptr
      top = bfi
      fc(7,bfi) = 0
      ptr = bf(j,-fc(5,bfi))

      do

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

        if ( ptr == bfi ) then
          exit
        end if

      end do

    end if
!
!  Find a boundary face visible from vertex I.
!
    do while ( top /= 0 )

      ptr = top
      top = fc(7,ptr)
      va = vm(fc(1,ptr))
      vb = vm(fc(2,ptr))
      vc = vm(fc(3,ptr))
      op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, vcl(1,vi) )

      if ( op == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS3 - Error!'
        write ( *, '(a)' ) '  Unexpected return value from OPSIDE.'
        ierr = 301
        return
      end if

      if ( op == 1 ) then

        front = ptr

        do while ( top /= 0 )

          ptr = top
          top = fc(7,ptr)
          fc(7,ptr) = -1

        end do

      else

        fc(7,ptr) = topnv
        topnv = ptr

      end if

    end do

    if ( front == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a)' ) '  FRONT = 0.'
      ierr = 306
      return
    end if
!
!  Find remaining visible boundary faces, add new tetrahedra with
!  vertex I, apply local transformation based on empty sphere criterion.
!
    call vbfac ( vcl(1,vi), ctr, vcl, vm, bf, fc, front, topnv )

    call nwthou ( i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, ntetra, &
      hdavbf, hdavfc, front, back, bfi, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a,i6)' ) '  NWTHOU IERR = ', ierr
      return
    end if

    call swapes ( .false., i, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
      ntetra, hdavfc, front, back, j, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a,i6)' ) '  SWAPES returned IERR = ', ierr
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

  return
end
