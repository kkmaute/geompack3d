subroutine bcdtri ( rflag, bf_num, nbpt, nipt, sizht, fc_max, vcl, vm, nfc, &
  nface, ntetra, fc, ht, ierr )

!*****************************************************************************80
!
!! BCDTRI constructs a boundary-constrained Delaunay triangulation in 3D.
!
!  Discussion:
!
!    This routine constructs the boundary-constrained Delaunay triangulation
!    of a set of 3D vertices, based on the local empty circumsphere criterion,
!    by using incremental approach and local transformations.
!
!    Vertices in interior of the convex hull are inserted one at a time in
!    the order given by the end of the VM array.  The initial tetrahedra
!    created due to a new vertex are obtained by a walk through the
!    triangulation until the location of the vertex is known.
!
!    If there are no interior vertices specified, then one may be added if
!    needed to produce a boundary-constrained triangulation.
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
!    Input, logical RFLAG, is TRUE iff return immediately (with input unchanged)
!    when NIPT = 0 and extra mesh vertex not removed.
!
!    Input, integer ( kind = 4 ) BF_NUM, the number of boundary faces or triangles.
!
!    Input, integer ( kind = 4 ) NBPT, the number of vertices (points) on boundary of
!    convex hull.
!
!    Input, integer ( kind = 4 ) NIPT, the number of vertices (points) in interior of
!    convex hull.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT; a good choice is a
!    prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) VM(1:NPT).  On input, indices of vertices of VCL
!    being triangulated where NPT = NBPT + MAX ( NIPT, 1 ); VM(1:NBPT)
!    are boundary points, rest are interior points; if NIPT = 0 then
!    VM(NPT) is an interior point which may be added to triangulation;
!    interior points are inserted in order VM(NBPT+1:NPT).
!    On output, VM(NPT) is set to 0 if NIPT = 0 and extra point not needed.
!
!    Input, integer ( kind = 4 ) FC(1:3,1:BF_NUM) - boundary triangles desired in triangulation;
!    entries are local vertex indices 1:NBPT (indices of VM)
!
!    Output, integer ( kind = 4 ) NFC, the number of positions used in FC array;
!    NFC <= FC_MAX.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation; NFACE <= NFC.
!
!    Output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:NFC), the array of face records which are in
!    linked lists in hash table with direct chaining. Fields are:
!    FC(1:3,*) - A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
!      vertices of face; if A <= 0, record is not used (it is
!      in linked list of avail records with indices <= NFC);
!      internal use: if B <= 0, face in queue, not in triang
!    FC(4:5,*) - D,E; indices in VM of 4th vertex of 1 or 2
!      tetrahedra with face ABC; if ABC is boundary face
!      then E = -1 (note that -E does not point to BF as in
!      routine DTRIW3 since array BF is not needed)
!    FC(6,*) - HTLINK; pointer (index in FC) of next element
!      in linked list (or NULL = 0)
!    FC(7,*) - used internally for QLINK (link for queues or
!      stacks); pointer (index in FC) of next face in queue/
!      stack (or NULL = 0); QLINK = -1 indicates face is not
!      in any queue/stack, and is output value (for records
!      not in avail list), except:
!    FC(7,2) - HDAVFC : head pointer of avail list in FC.
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and tetrahedra of triangulation.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) nbpt
  integer ( kind = 4 ) nipt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(1)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) ptr
  logical              remov
  logical              rflag
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(nbpt+nipt)
!
!  Initialize triangulation as first interior vertex joined to all
!  boundary faces.
!
  ierr = 0

  if ( fc_max < 5 * bf_num / 2 ) then
    ierr = 11
    return
  end if

  ht(0:sizht-1) = 0

  npt = nbpt + max ( nipt, 1 )
  d = nbpt + 1
  hdavfc = 0
  nfc = bf_num
  ntetra = bf_num

  do i = 1, bf_num
    call htins(i,fc(1,i),fc(2,i),fc(3,i),d,-1,npt,sizht,fc,ht)
  end do

  if ( msglvl == 4 ) then
    write ( *,600) bf_num,nbpt,nipt
  end if

  do i = 1, bf_num

    a = fc(1,i)
    b = fc(2,i)
    c = fc(3,i)

    if ( msglvl == 4 ) then
      write ( *,610) a,b,c,d,vm(a),vm(b),vm(c), vm(d)
    end if

    do j = 1, 3

      if ( j == 2 ) then
        b = fc(3,i)
        c = fc(2,i)
      else if ( j == 3 ) then
        a = fc(2,i)
        c = fc(1,i)
      end if

      ind = htsrc ( a, b, d, npt, sizht, fc, ht )

      if ( ind <= 0 ) then
        nfc = nfc + 1
        call htins(nfc,a,b,d,c,0,npt,sizht,fc,ht)
      else
        fc(5,ind) = c
      end if

    end do

  end do
!
!  If NIPT=0, apply local transformations to try to remove extra point.
!  Then apply local transformations based on local empty sphere criterion.
!  This latter step is also performed for 0 < NIPT.  Note that BF is
!  a dummy array in call to SWAPES since BF is not referenced.
!
  if ( nipt == 0 ) then

    call swprem ( bf_num, nbpt, sizht, fc_max, nfc, ntetra, vcl, vm, fc, &
      ht, remov, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( remov ) then
      vm(npt) = 0
    end if

    if ( .not. remov .and. rflag ) then
      return
    end if
    hdavfc = fc(7,2)
    fc(7,2) = -1
  end if

  front = 0

  do i = bf_num+1, nfc

    if ( 0 < fc(1,i) ) then
      if ( front == 0 ) then
        front = i
      else
        fc(7,back) = i
      end if
      back = i
    end if

  end do

  if ( front /= 0 ) then
    fc(7,back) = 0
  end if

  if ( msglvl == 4 ) then
    write ( *,620) d,vm(d)
  end if

  call swapes ( .true., 0, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
    ntetra, hdavfc, front, back, ind, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( nipt <= 1 ) then
    go to 90
  end if
!
!  Insert I-th vertex into pseudo-locally optimal triangulation of first I-1
!  vertices.  Walk through triangulation to find location of vertex I, create
!  new tetrahedra, apply local transformations based on empty sphere criterion.
!
  a = 1
  b = 0
  ifac = nfc

  do while ( fc(1,ifac) <= 0 )
    ifac = ifac - 1
  end do

  do i = nbpt+2, npt

    vi = vm(i)

    if ( msglvl == 4 ) then
      write ( *,620) i, vi
    end if

    if ( fc(5,ifac) == i - 1 ) then
      ivrt = 5
    else
      ivrt = 4
    end if

    call walkt3(vcl(1,vi),npt,sizht,ntetra,vcl,vm,fc,ht,ifac,ivrt, ierr )

    if ( ierr == 307 ) then

      ierr = 0

      if ( 2 <= msglvl ) then
        write ( *,630) i
      end if

      call lsrct3(vcl(1,vi),npt,sizht,nfc,vcl,vm,fc,ht,ifac,ivrt, ierr )

      if ( ifac == 0 ) then
        ierr = 331
      end if

    end if

    if ( ierr /= 0 ) then
      return
    end if

    if ( ivrt == 6 ) then

      if ( fc(5,ifac) <= 0 ) then
        ierr = 331
      else
        call nwthfc(i,ifac,npt,sizht,a,nfc,1,fc_max,bf,fc,ht, &
          ntetra,b,hdavfc,front,back, ierr )
      end if

    else if ( 4 <= ivrt ) then

      call nwthin(i,ifac,ivrt,npt,sizht,nfc,fc_max,fc,ht,ntetra, &
        hdavfc,front,back, ierr )

    else if ( 1 <= ivrt ) then

      call nwthed(i,ifac,ivrt,npt,sizht,a,nfc,1,fc_max,bf,fc,ht, &
        ntetra,b,hdavfc,front,back, ierr )

      if ( ierr == 12 ) then
        ierr = 331
      end if

    else

      ierr = 331

    end if

    if ( ierr /= 0 ) then
      return
    end if

    call swapes ( .true., i, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
      ntetra, hdavfc, front, back, ind, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( ind /= 0 ) then
      ifac = ind
    end if

  end do
!
!  Make final pass based on local empty circumsphere criterion.
!
  front = 0

  do i = bf_num+1, nfc
    if ( 0 < fc(1,i) ) then
      if ( front == 0 ) then
        front = i
      else
        fc(7,back) = i
      end if
      back = i
    end if
  end do

  if ( front /= 0 ) then
    fc(7,back) = 0
  end if

  call swapes ( .true., 0, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
    ntetra, hdavfc, front, back, ind, ierr )

  if ( ierr /= 0 ) then
    return
  end if

   90 continue

  nface = nfc
  ptr = hdavfc

  do while ( ptr /= 0 )
    nface = nface - 1
    ptr = -fc(1,ptr)
  end do

  fc(7,2) = hdavfc

  600 format (/' bcdtri: bf_num,nbpt,nipt=',3i7)
  610 format (1x,4i7,' : ',4i7)
  620 format (/1x,'step',i7,':   vertex i =',i7)
  630 format (1x,'linear search required in step',i7)

  return
end
