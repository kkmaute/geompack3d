subroutine dtriw3 ( npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, nface, &
  ntetra, bf, fc, ht, ierr )

!*****************************************************************************80
!
!! DTRIW3 constructs a Delaunay triangulation of vertices in 3D.
!
!  Discussion:
!
!    This routine constructs a Delaunay triangulation of 3D vertices using
!    incremental approach and local transformations.  Vertices are
!    inserted one at a time in order given by VM array.  The initial
!    tetrahedra created due to a new vertex are obtained by a walk
!    through the triangulation until location of vertex is known.
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
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT; a good choice is a
!    prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input/output, VM(1:NPT), the indices of vertices of VCL being
!    triangulated; vertices are inserted in order given by VM, except that
!    on output, the third and fourth elements may be swapped so that first
!    4 vertices are non-coplanar.
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
!    Output, integer ( kind = 4 ) BF(1:3,1:BF_NUM), the array of boundary face records
!    containing pointers (indices) to FC; if FC(5,I) = -J < 0 and
!    FC(1:3,I) = ABC, then BF(1,J) points to other boundary face with edge BC,
!    BF(2,J) points to other boundary face with edge AC, and
!    BF(3,J) points to other boundary face with edge AB;
!    if BF(1,J) <= 0, record is not used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:NFC), the array of face records which
!    are in linked lists in hash table with direct chaining. Fields are:
!    FC(1:3,*) - A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
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
!    FC(7,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC.
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and tetrahedra of triangulation
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  real    ( kind = 8 ) ctr(3)
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)

  ierr = 0
!
!  Reorder points so that first four points are in general position.
!
  call frstet ( .false., npt, vcl, vm, i3, i4, ierr )

  if ( ierr /= 0 ) then
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
  call htins(1,1,2,3,4,-1,npt,sizht,fc,ht)
  call htins(2,1,2,4,3,-2,npt,sizht,fc,ht)
  call htins(3,1,3,4,2,-3,npt,sizht,fc,ht)
  call htins(4,2,3,4,1,-4,npt,sizht,fc,ht)
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
  ifac = 4

  if ( msglvl == 4 ) then
    write ( *,600) (vm(i),i=1,4),i3,i4
  end if
!
!  Insert Ith vertex into Delaunay triang of first I-1 vertices.
!  Walk through triang to find location of vertex I, create new
!  tetrahedra, apply local transf based on empty sphere criterion.
!
  do i = 5, npt

    vi = vm(i)

    if ( msglvl == 4 ) then
      write ( *,610) i,vi
    end if

    if ( fc(5,ifac) == i-1 ) then
      ivrt = 5
    else
      ivrt = 4
    end if

    call walkt3(vcl(1,vi),npt,sizht,ntetra,vcl,vm,fc,ht,ifac,ivrt, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( ivrt == 6 ) then

      call nwthfc(i,ifac,npt,sizht,bf_num,nfc,bf_max,fc_max,bf,fc,ht, &
        ntetra,hdavbf,hdavfc,front,back, ierr )

    else if ( 4 <= ivrt ) then

      call nwthin(i,ifac,ivrt,npt,sizht,nfc,fc_max,fc,ht,ntetra, &
        hdavfc,front,back, ierr )

    else if ( ivrt == 0 ) then

      front = ifac
      call vbfac ( vcl(1,vi), ctr, vcl, vm, bf, fc, front, 0 )

      if ( ierr /= 0 ) then
        return
      end if

      call nwthou ( i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, ntetra, &
        hdavbf, hdavfc, front, back, ind, ierr )

    else if ( 1 <= ivrt ) then

      call nwthed(i,ifac,ivrt,npt,sizht,bf_num,nfc,bf_max,fc_max,bf,fc, &
        ht,ntetra,hdavbf,hdavfc,front,back, ierr )

    else

      ierr = 302

    end if

    if ( ierr /= 0 ) then
      return
    end if

    call swapes ( .false., i, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
      ntetra, hdavfc, front, back, ind, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( ind /= 0 ) then
      ifac = ind
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

  600 format (/1x,'dtriw3: first tetrahedron: ',4i7/4x,'i3, i4 =',2i7)
  610 format (/1x,'step',i7,':   vertex i =',i7)

  return
end
