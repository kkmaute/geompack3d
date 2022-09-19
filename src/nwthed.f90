subroutine nwthed ( i, ifac, iedg, npt, sizht, bf_num, nfc, bf_max, fc_max, &
  bf, fc, ht, ntetra, hdavbf, hdavfc, front, back, ierr )

!*****************************************************************************80
!
!! NWTHED creates new tetrahedra from the insertion of a vertex, in 3D.
!
!  Discussion:
!
!    This routine creates new tetrahedra in a 3D triangulation from the
!    insertion of vertex I on edge FC(IEDG,IFAC).
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
!    Input, integer ( kind = 4 ) I, the (local) index of next vertex inserted in
!    triangulation; it is assumed I is largest index so far.
!
!    Input, integer ( kind = 4 ) IFAC, IEDG, the edge containing I is FC(IEDG,IFAC);
!    1 <= IEDG <= 3.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
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
!    Input/output, integer ( kind = 4 ) BF(1:3,1:BF_MAX), the array of boundary face records;
!    see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVBF, the head pointer to available BF records.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of queue
!    of interior faces ABC such that ABCI is a new tetrahedron.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bb
  logical              bedge
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfn(2)
  integer ( kind = 4 ) bfo(2)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cp
  integer ( kind = 4 ) csav
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fcn(2)
  integer ( kind = 4 ) fco(2)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) inew
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbrac
  integer ( kind = 4 ) nbrbc
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) nfcin
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra

  ierr = 0
  front = 0
  back = 0

  if ( iedg == 1 ) then
    a = fc(1,ifac)
    b = fc(2,ifac)
    c = fc(3,ifac)
  else if ( iedg == 2 ) then
    a = fc(2,ifac)
    b = fc(3,ifac)
    c = fc(1,ifac)
  else
    a = fc(1,ifac)
    b = fc(3,ifac)
    c = fc(2,ifac)
  end if

  csav = c
  d = fc(4,ifac)
!
!  Determine faces incident on edge AB in circular order, and store
!  their indices at end of FC array.
!
  bedge = .false.
  ind = ifac
  k = 0

  do

    k = k + 1

    if ( fc_max < nfc + k ) then
      ierr = 11
      return
    end if

    fc(1,nfc+k) = ind

    if ( bedge ) then
      exit
    end if

    if ( d == csav ) then
      go to 50
    end if

    ind = htsrc(a,b,d,npt,sizht,fc,ht)

    if ( ind <= 0 ) then
      ierr = 300
      return
    end if

    if ( fc(5,ind) <= 0 ) then
      bedge = .true.
      cp = d
    else if ( fc(4,ind) == c ) then
      c = d
      d = fc(5,ind)
    else
      c = d
      d = fc(4,ind)
    end if

  end do
!
!  Determine further faces in case of AB being a boundary edge.
!
20 continue

  if ( fc(5,ifac) <= 0 ) then
    go to 50
  end if

  l = k

  do j = 1, k/2
    e = fc(1,nfc+j)
    fc(1,nfc+j) = fc(1,nfc+l)
    fc(1,nfc+l) = e
    l = l - 1
  end do

  c = csav
  csav = cp
  d = fc(5,ifac)

  do

    ind = htsrc(a,b,d,npt,sizht,fc,ht)

    if ( ind <= 0 ) then
      ierr = 300
      return
    end if

    k = k + 1

    if ( fc_max < nfc + k ) then
      ierr = 11
      return
    end if

    fc(1,nfc+k) = ind

    if ( fc(5,ind) <= 0 ) then
      exit
    else if ( fc(4,ind) == c ) then
      c = d
      d = fc(5,ind)
    else
      c = d
      d = fc(4,ind)
    end if

  end do
!
!  Create new faces and tetrahedra, and add faces to queue.
!
50 continue

  nfcin = nfc
  nfc = nfc + k
  ntetra = ntetra + k

  if ( bedge ) then
    ntetra = ntetra - 1
    fcn(1) = nfcin + 1
    fcn(2) = nfcin + k
    fco(1) = fc(1,nfcin+1)
    fco(2) = fc(1,nfcin+k)
    bfo(1) = -fc(5,fco(1))
    bfo(2) = -fc(5,fco(2))
  end if

  do j = 1, k

    inew = nfcin + j
    ind = fc(1,inew)

    if ( fc(1,ind) == a ) then
      if ( fc(2,ind) == b ) then
        c = fc(3,ind)
      else
        c = fc(2,ind)
      end if
    else
      c = fc(1,ind)
    end if

    d = fc(4,ind)
    e = fc(5,ind)
    call htdel(ind,npt,sizht,fc,ht)
    call htins(ind,a,c,i,d,e,npt,sizht,fc,ht)
    call htins(inew,b,c,i,d,e,npt,sizht,fc,ht)

    if ( j == k ) then
      if ( bedge ) then
        cycle
      end if
      d = csav
    else if ( 1 < j ) then
      if ( d == cp ) then
        d = e
      end if
    end if

    call availf ( hdavfc, nfc, fc_max, fc, ind, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call htins(ind,c,d,i,a,b,npt,sizht,fc,ht)
    aa = a
    bb = b

    do l = 1, 2

      if ( l == 2 ) then
        aa = b
        bb = a
      end if

      ind = htsrc(aa,c,d,npt,sizht,fc,ht)

      if ( ind <= 0 ) then
        ierr = 300
        return
      end if

      if ( fc(5,ind) <= 0 ) then
        fc(4,ind) = i
      else
        if ( fc(4,ind) == bb ) then
          fc(4,ind) = i
        else
          fc(5,ind) = i
        end if
        if ( front == 0 ) then
          front = ind
        else
          fc(7,back) = ind
        end if
        back = ind
      end if

      if ( msglvl == 4 ) then
        write ( *,600) aa,c,d,i
      end if

    end do

    cp = c

  end do

  if ( front /= 0 ) then
    fc(7,back) = 0
  end if

  if ( bedge ) then

    d = c
    c = csav

    if ( hdavbf /= 0 ) then

      bfn(1) = hdavbf
      hdavbf = -bf(1,hdavbf)

      if ( hdavbf /= 0 ) then
        bfn(2) = hdavbf
        hdavbf = -bf(1,hdavbf)
      else
        bf_num = bf_num + 1
        bfn(2) = bf_num
      end if

    else

      bf_num = bf_num + 2
      bfn(1) = bf_num - 1
      bfn(2) = bf_num

    end if

    if ( bf_max < bf_num ) then
      bf_num = bf_max
      ierr = 12
      return
    end if

    fc(5,nfcin+1) = -bfn(1)
    fc(5,nfcin+k) = -bfn(2)

    do j = 1, 2

      if ( j == 2 ) then
        c = d
      end if

      if ( c < a ) then
        nbrac = bf(3,bfo(j))
        nbrbc = bf(2,bfo(j))
        bf(1,bfo(j)) = fco(3-j)
        bf(2,bfo(j)) = fcn(j)
        bf(1,bfn(j)) = fcn(3-j)
        bf(2,bfn(j)) = fco(j)
      else if ( c < b ) then
        nbrac = bf(3,bfo(j))
        nbrbc = bf(1,bfo(j))
        bf(1,bfo(j)) = fcn(j)
        bf(2,bfo(j)) = fco(3-j)
        bf(1,bfn(j)) = fcn(3-j)
        bf(2,bfn(j)) = fco(j)
      else
        nbrac = bf(2,bfo(j))
        nbrbc = bf(1,bfo(j))
        bf(1,bfo(j)) = fcn(j)
        bf(2,bfo(j)) = fco(3-j)
        bf(1,bfn(j)) = fco(j)
        bf(2,bfn(j)) = fcn(3-j)
      end if

      bf(3,bfo(j)) = nbrac
      bf(3,bfn(j)) = nbrbc
      l = -fc(5,nbrbc)

      if ( bf(1,l) == fco(j) ) then
        bf(1,l) = fcn(j)
      else if ( bf(2,l) == fco(j) ) then
        bf(2,l) = fcn(j)
      else
        bf(3,l) = fcn(j)
      end if

    end do

  end if

  600 format ( '  New tetra: ',4i7)

  return
end
