subroutine swprem ( bf_num, nbpt, sizht, fc_max, nfc, ntetra, vcl, vm, fc, ht, &
  remov, ierr )

!*****************************************************************************80
!
!! SWPREM tries to remove an interior vertex.
!
!  Discussion:
!
!    This routine tries to remove an interior vertex in a boundary-constrained
!    3D triangulation in which the interior vertex is joined to all boundary
!    faces.
!
!    Local transformations are applied to try to remove the interior vertex,
!    i.e. get a boundary-constrained triangulation formed from the given
!    boundary faces only.
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
!    Input, integer ( kind = 4 ) BF_NUM, the number of boundary faces or triangles.
!
!    Input, integer ( kind = 4 ) NBPT, the number of vertices (points) on boundary of
!    convex hull.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array;
!    NFC <= FC_MAX.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NBPT+1), the indices of vertices of VCL being
!    triangulated where first NBPT are boundary points and last is interior.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:NFC), the array of face records; see
!    routine BCDTRI; first NBPT columns are for boundary faces, rest are
!    interior.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) FC(7,2), the HDAVFC : head pointer of avail list in FC.
!
!    Output, logical REMOV, TRUE iff interior point VM(NBPT+1) is removed.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  real    ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  logical              degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) itet
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbpt
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) npt
  logical              remov
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(nbpt+1)
!
!  Put all interior faces in queue for processing. I is index of
!  interior point. ITET is number of tetrahedra incident on I.
!
  ierr = 0

  if ( msglvl == 4 ) then
    write ( *,650)
  end if

  npt = nbpt + 1
  i = npt
  vi = vm(i)
  hdavfc = 0
  itet = bf_num

  if ( itet == 4 ) then
    go to 40
  end if

  front = bf_num + 1

  do j = bf_num+1, nfc-1
    fc(7,j) = j + 1
  end do

  fc(7,nfc) = 0
  back = nfc

20 continue

  if ( front == 0 ) then
    go to 40
  end if

  ind = front
  front = fc(7,ind)

  if ( fc(2,ind) == 0 ) then

    if ( ind == nfc ) then
      nfc = nfc - 1
    else
      fc(1,ind) = -hdavfc
      hdavfc = ind
    end if

    go to 20

  end if

  fc(7,ind) = -1
  a = fc(1,ind)
  b = fc(2,ind)
  c = fc(4,ind)
  d = fc(5,ind)
  va = vm(a)
  vb = vm(b)
  vc = vm(c)
  vd = vm(d)

  if ( msglvl == 4 ) then
    write ( *,600) ind,a,b,i,c,d
  end if

  call baryth ( vcl(1,va), vcl(1,vb), vcl(1,vi), vcl(1,vc), vcl(1,vd), &
    alpha, degen )

  if ( degen ) then
    ierr = 301
    return
  else if ( 0.0D+00 < alpha(4) ) then
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
!
!  Swap 2 tetrahedra for 3.
!
  if ( kneg == 1 .and. kzero == 0 ) then

    call updatr ( a, b, c, i, d, .false., npt, sizht, front, back, fc, &
      ht, ierr )

    call updatr ( a, b, d, i, c, .false., npt, sizht, front, back, fc, &
      ht, ierr )

    call updatr(a,c,i,b,d,.true.,npt,sizht,front,back,fc,ht, ierr )
    call updatr(b,c,i,a,d,.true.,npt,sizht,front,back,fc,ht, ierr )
    call updatr(a,d,i,b,c,.true.,npt,sizht,front,back,fc,ht, ierr )
    call updatr(b,d,i,a,c,.true.,npt,sizht,front,back,fc,ht, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call htdel(ind,npt,sizht,fc,ht)
    call htins(ind,a,c,d,b,i,npt,sizht,fc,ht)
    call availf(hdavfc,nfc,fc_max,fc,ind,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    call htins(ind,b,c,d,a,i,npt,sizht,fc,ht)
    call availf(hdavfc,nfc,fc_max,fc,ind,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    call htins(ind,c,d,i,a,b,npt,sizht,fc,ht)
    ntetra = ntetra + 1

    if ( msglvl == 4 ) then
      write ( *,610)
    end if
!
!  Relabel so that AI intersects BCD.
!
  else if ( kneg == 2 .and. kzero == 0 ) then

    if ( alpha(3) < 0.0D+00 ) go to 20

    if ( alpha(1) < 0.0D+00 ) then
      call i4_swap ( a, b )
    end if

    ind1 = htsrc(a,c,i,npt,sizht,fc,ht)

    if ( ind1 <= 0 ) then
      ierr = 300
      return
    end if
!
!  Swap 3 tetrahedra for 2.
!
    if ( fc(4,ind1) == d .or. fc(5,ind1) == d ) then

      va = vm(a)
      vb = vm(b)
      call updatr(a,b,c,i,d,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(a,b,d,i,c,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(a,c,d,i,b,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(b,c,i,a,d,.true.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(b,d,i,a,c,.true.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(c,d,i,a,b,.true.,npt,sizht,front,back,fc,ht, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htdel(ind,npt,sizht,fc,ht)
      call htins(ind,b,c,d,a,i,npt,sizht,fc,ht)
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

      ind1 = htsrc(a,d,i,npt,sizht,fc,ht)

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
      itet = itet - 2

      if ( msglvl == 4 ) then
        write ( *,620) b,c,d
      end if

    end if

  else if ( kneg == 1 .and. kzero == 1 ) then

    if ( alpha(3) == 0.0D+00) go to 20
!
!  Relabel vertices so that BI intersects CD.
!
    if ( alpha(2) == 0.0D+00 ) then
      call i4_swap ( a, b )
    end if

    ind1 = htsrc(b,c,i,npt,sizht,fc,ht)
    ind2 = htsrc(b,d,i,npt,sizht,fc,ht)

    if ( ind1 <= 0 .or. ind2 <= 0 ) then
      ierr = 300
      return
    end if

    if ( fc(4,ind1) == a ) then
      e = fc(5,ind1)
    else
      e = fc(4,ind1)
    end if

    if ( fc(4,ind2) == a ) then
      f = fc(5,ind2)
    else
      f = fc(4,ind2)
    end if
!
!  Swap 4 tetrahedra for 4.
!
    if ( e == f ) then

      va = vm(a)
      vb = vm(b)
      ve = vm(e)
      call updatr(a,b,c,i,d,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(a,b,d,i,c,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(b,c,e,i,d,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(b,d,e,i,c,.false.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(a,c,i,b,d,.true.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(a,d,i,b,c,.true.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(c,e,i,b,d,.true.,npt,sizht,front,back,fc,ht, ierr )
      call updatr(d,e,i,b,c,.true.,npt,sizht,front,back,fc,ht, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htdel(ind,npt,sizht,fc,ht)
      call htins(ind,a,c,d,b,i,npt,sizht,fc,ht)
      ind = htsrc(b,e,i,npt,sizht,fc,ht)

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

      call htins(ind,c,d,e,b,i,npt,sizht,fc,ht)
      call htdel(ind1,npt,sizht,fc,ht)

      if ( 0 <= fc(7,ind1) ) then
        fc(2,ind1) = 0
        call availf(hdavfc,nfc,fc_max,fc,ind1,ierr)

        if ( ierr /= 0 ) then
          return
        end if

      end if

      call htins(ind1,b,c,d,a,e,npt,sizht,fc,ht)
      call htdel(ind2,npt,sizht,fc,ht)
      j = fc(7,ind2)
      call htins(ind2,c,d,i,a,e,npt,sizht,fc,ht)
      fc(7,ind2) = j

      if ( fc(7,ind2) == -1 ) then

        fc(7,ind2) = 0

        if ( front == 0 ) then
          front = ind2
        else
          fc(7,back) = ind2
        end if

        back = ind2

      end if

      itet = itet - 2

      if ( msglvl == 4 ) then
        write ( *,630) b,i,c,d,e
      end if

    end if

  end if

  go to 20

40 continue

  if ( 2 <= msglvl ) then
    write ( *,660) itet,bf_num,nbpt
  end if

  if ( itet /= 4 ) then
    remov = .false.
  else
!
!  Replace tetrahedra ABCI, ABDI, ACDI, BCDI by ABCD.
!
    remov = .true.
    ind = nbpt + 1

50  continue

    if ( fc(3,ind) /= i .or. fc(1,ind) <= 0 ) then
      ind = ind + 1
      go to 50
    end if

    a = fc(1,ind)
    b = fc(2,ind)
    c = fc(4,ind)
    d = fc(5,ind)
    call updatr(a,b,c,i,d,.false.,npt,sizht,front,back,fc,ht, ierr )
    call updatr(a,b,d,i,c,.false.,npt,sizht,front,back,fc,ht, ierr )
    call updatr(a,c,d,i,b,.false.,npt,sizht,front,back,fc,ht, ierr )
    call updatr(b,c,d,i,a,.false.,npt,sizht,front,back,fc,ht, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    aa = a
    bb = c

    do j = 1,6

      call htdel(ind,npt,sizht,fc,ht)

      if ( ind == nfc ) then
        nfc = nfc - 1
      else
        fc(1,ind) = -hdavfc
        hdavfc = ind
      end if

      if ( j == 6 ) then
        cycle
      end if

      if ( j == 2 ) then
        bb = d
      else if ( j == 3 ) then
        aa = b
      else if ( j == 4 ) then
        bb = c
      else if ( j == 5 ) then
        aa = d
      end if

      ind = htsrc ( aa, bb, i, npt, sizht, fc, ht )

      if ( ind <= 0 ) then
        ierr = 300
        return
      end if

    end do

    ntetra = ntetra - 3

    if ( msglvl == 4 ) then
      write ( *,640) i,a,b,c,d
    end if

  end if

  fc(7,2) = hdavfc

  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   e =',i7)
  640 format (1x,'delete ',i7,' from tetra ',4i7)
  650 format (/1x,'swprem')
  660 format (1x,'swprem: itet,bf_num,nbpt=',3i5)

  return
end
