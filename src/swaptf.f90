subroutine swaptf ( top, npt, sizht, nfc, fc_max, vcl, vm, fc, ht, ntetra, &
  hdavfc, front, back, ierr )

!*****************************************************************************80
!
!! SWAPTF swaps tranformable faces in a 3D triangulation.
!
!  Discussion:
!
!    This routine swaps transformable faces of type T23, T32, or T44 in a 3D
!    triangulation in the order given on stack.
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
!    Input/output, integer ( kind = 4 ) TOP, the pointer to stack of faces to be
!    transformed; each face should be transformable upon reaching top of
!    stack.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of
!    queue of interior faces for which local optimality test may be applied;
!    this routine adds faces to end of queue.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  real    ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  logical              degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) top
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vm(npt)

  ierr = 0

10 continue

  if ( top == 0 ) then
    return
  end if

  ind = top
  top = fc(7,ind)

  if ( fc(2,ind) == 0 ) then
    ierr = 308
    return
  end if

  fc(7,ind) = -1
  a = fc(1,ind)
  b = fc(2,ind)
  c = fc(3,ind)
  d = fc(4,ind)
  e = fc(5,ind)
  va = vm(a)
  vb = vm(b)
  vc = vm(c)
  vd = vm(d)
  ve = vm(e)

  if ( msglvl == 4 ) then
    write ( *,600) ind,a,b,c,d,e
  end if

  call baryth(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve), &
    alpha,degen)

  if ( degen ) then
    ierr = 301
    return
  end if

  if ( 0.0D+00 < alpha(4) ) then
    ierr = 309
    return
  end if

  kneg = 1
  kzero = 0

  do j = 1, 3
    if ( alpha(j) < 0.0D+00 ) then
      kneg = kneg + 1
    else if ( alpha(j) == 0.0D+00 ) then
      kzero = kzero + 1
    end if
  end do

  if ( kneg == 1 .and. kzero == 0 ) then
!
!  Swap 2 tetrahedra for 3.
!
    call updatf(a,b,d,c,e,0,npt,sizht,front,back,fc,ht, ierr )
    call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
    call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
    call updatf(a,b,e,c,d,0,npt,sizht,front,back,fc,ht, ierr )
    call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
    call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call htdel(ind,npt,sizht,fc,ht)
    call htins(ind,a,d,e,b,c,npt,sizht,fc,ht)
    call availf(hdavfc,nfc,fc_max,fc,ind,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    call htins(ind,b,d,e,a,c,npt,sizht,fc,ht)
    call availf(hdavfc,nfc,fc_max,fc,ind,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
    ntetra = ntetra + 1

    if ( msglvl == 4 ) then
      write ( *,610)
    end if

  else if ( kneg == 2 .and. kzero == 0 ) then
!
!  Relabel so new tetrahedra are ACDE, BCDE (edge AB deleted).
!
    if ( alpha(1) < 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) < 0.0D+00 ) then
      call i4_swap ( b, c )
    end if

    ind1 = htsrc(a,b,d,npt,sizht,fc,ht)

    if ( ind1 <= 0 ) then
      ierr = 300
      return
    end if
!
!  Swap 3 tetrahedra for 2.
!
    if ( fc(4,ind1) == e .or. fc(5,ind1) == e ) then

      call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,d,e,b,c,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,d,e,a,c,0,npt,sizht,front,back,fc,ht, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htdel(ind,npt,sizht,fc,ht)
      call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
      call htdel(ind1,npt,sizht,fc,ht)

      if ( 0 <= fc(7,ind1) ) then
        fc(2,ind1) = 0
      else
        if ( ind1 == nfc ) then
          nfc = nfc - 1
        else
          fc(1,ind1) = -hdavfc
          hdavfc = ind1
        end if
      end if

      ind1 = htsrc(a,b,e,npt,sizht,fc,ht)

      if ( ind1 <= 0 ) then
        ierr = 300
        return
      end if

      call htdel(ind1,npt,sizht,fc,ht)

      if ( 0 <= fc(7,ind1) ) then
        fc(2,ind1) = 0
      else
        if ( ind1 == nfc ) then
          nfc = nfc - 1
        else
          fc(1,ind1) = -hdavfc
          hdavfc = ind1
        end if
      end if

      ntetra = ntetra - 1

      if ( msglvl == 4 ) then
        write ( *,620) c,d,e
      end if

    else

      ierr = 308
      return

    end if

  else if ( kneg == 1 .and. kzero == 1 ) then
!
!  Relabel vertices so that DE intersects AB.
!
    if ( alpha(1) == 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) == 0.0D+00 ) then
      call i4_swap ( b, c )
    end if

    ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
    ind2 = htsrc(a,b,e,npt,sizht,fc,ht)

    if ( ind1 <= 0 .or. ind2 <= 0 ) then
      ierr = 300
      return
    end if

    if ( fc(4,ind1) == c ) then
      f = fc(5,ind1)
    else
      f = fc(4,ind1)
    end if

    if ( fc(4,ind2) == c ) then
      g = fc(5,ind2)
    else
      g = fc(4,ind2)
    end if

    if ( f <= 0 .or. g <= 0 .or. f /= g ) then
      ierr = 308
      return
    else
!
!  Swap 4 tetrahedra for 4.
!
      call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,d,f,b,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,e,f,b,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,d,f,a,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,e,f,a,d,0,npt,sizht,front,back,fc,ht, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htdel(ind,npt,sizht,fc,ht)
      call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
      ind = htsrc(a,b,f,npt,sizht,fc,ht)

      if ( ind <= 0 ) then
        ierr = 300
        return
      end if

      call htdel(ind,npt,sizht,fc,ht)

      if ( 0 <= fc(7,ind) ) then
        fc(2,ind) = 0
        call availf(hdavfc,nfc,fc_max,fc,ind,ierr)
        if ( ierr /= 0 ) then
          return
        end if
      end if

      call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)

      call htdel(ind1,npt,sizht,fc,ht)

      if ( 0 <= fc(7,ind1) ) then

        fc(2,ind1) = 0
        call availf(hdavfc,nfc,fc_max,fc,ind1,ierr)

        if ( ierr /= 0 ) then
          return
        end if

      end if

      call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)

      call htdel(ind2,npt,sizht,fc,ht)

      if ( 0 <= fc(7,ind2) ) then
        fc(2,ind2) = 0
        call availf(hdavfc,nfc,fc_max,fc,ind2,ierr)
        if ( ierr /= 0 ) then
          return
        end if
      end if

      call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)

      if ( front == 0 ) then
        front = ind1
      else
        fc(7,back) = ind1
      end if

      fc(7,ind1) = ind2
      fc(7,ind2) = 0
      back = ind2

      if ( msglvl == 4 ) then
        write ( *,630) a,b,d,e,f
      end if

    end if

  else

    ierr = 308
    return

  end if

  go to 10

  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)

end
