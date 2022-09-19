subroutine nwthfc ( i, ifac, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, &
  fc, ht, ntetra, hdavbf, hdavfc, front, back, ierr )

!*****************************************************************************80
!
!! NWTHFC creates new tetrahedra after the insertion of a new face vertex.
!
!  Discussion:
!
!    This routine creates new tetrahedra in a 3D triangulation from the
!    insertion of vertex I on face FC(*,IFAC).
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
!    Input, integer ( kind = 4 ) IFAC, the face containing vertex I is FC(*,IFAC).
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
!    Input/output, integer ( kind = 4 ) BF(1:3,1:BF_MAX), the array of boundary face
!    records; see DTRIS3.
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
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bf1
  integer ( kind = 4 ) bf2
  logical              bface
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbrac
  integer ( kind = 4 ) nbrbc
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) nv

  ierr = 0
  front = 0
  back = 0
  nv = 5
  bface = ( fc(5,ifac) <= 0 )

  if ( bface ) then
    nv = 4
  end if

  a = fc(1,ifac)
  b = fc(2,ifac)
  c = fc(3,ifac)

  do iv = 4, nv

    ntetra = ntetra + 2
    d = fc(iv,ifac)

    do j = 1, 3

      if ( j == 1 ) then
        aa = a
        bb = b
        cc = c
      else if ( j == 2 ) then
        aa = b
        bb = c
        cc = a
      else
        aa = c
        bb = a
        cc = b
      end if

      call availf ( hdavfc, nfc, fc_max, fc, ind, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htins(ind,aa,d,i,bb,cc,npt,sizht,fc,ht)
      ind = htsrc(aa,bb,d,npt,sizht,fc,ht)

      if ( ind <= 0 ) then
        ierr = 300
        return
      end if

      if ( fc(4,ind) == cc ) then
        fc(4,ind) = i
      else
        fc(5,ind) = i
      end if

      if ( 0 < fc(5,ind) ) then
        if ( front == 0 ) then
          front = ind
        else
          fc(7,back) = ind
        end if
        back = ind
      end if

      if ( msglvl == 4 ) then
        write ( *,600) aa,bb,d,i
      end if

    end do

  end do

  if ( front /= 0 ) then
    fc(7,back) = 0
  end if

  call availf ( hdavfc, nfc, fc_max, fc, ind, ierr )
  call availf ( hdavfc, nfc, fc_max, fc, ind2, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( bface ) then
    e = fc(5,ifac)
  else
    e = fc(4,ifac)
  end if

  call htdel(ifac,npt,sizht,fc,ht)
  call htins(ifac,a,b,i,d,e,npt,sizht,fc,ht)
  call htins(ind,a,c,i,d,e,npt,sizht,fc,ht)
  call htins(ind2,b,c,i,d,e,npt,sizht,fc,ht)

  if ( bface ) then

    e = -e

    if ( hdavbf /= 0 ) then
      bf1 = hdavbf
      hdavbf = -bf(1,hdavbf)
      if ( hdavbf /= 0 ) then
        bf2 = hdavbf
        hdavbf = -bf(1,hdavbf)
      else
        bf_num = bf_num + 1
        bf2 = bf_num
      end if
    else
      bf_num = bf_num + 2
      bf1 = bf_num - 1
      bf2 = bf_num
    end if

    if ( bf_max < bf_num ) then
      bf_num = bf_max
      ierr = 12
      return
    end if

    fc(5,ind) = -bf1
    fc(5,ind2) = -bf2
    nbrac = bf(2,e)
    nbrbc = bf(1,e)
    bf(1,e) = ind2
    bf(2,e) = ind
    bf(1,bf1) = ind2
    bf(2,bf1) = ifac
    bf(3,bf1) = nbrac
    bf(1,bf2) = ind
    bf(2,bf2) = ifac
    bf(3,bf2) = nbrbc
    j = -fc(5,nbrac)

    if ( bf(1,j) == ifac ) then
      bf(1,j) = ind
    else if ( bf(2,j) == ifac ) then
      bf(2,j) = ind
    else
      bf(3,j) = ind
    end if

    j = -fc(5,nbrbc)

    if ( bf(1,j) == ifac ) then
      bf(1,j) = ind2
    else if ( bf(2,j) == ifac ) then
      bf(2,j) = ind2
    else
      bf(3,j) = ind2
    end if

  end if

  600 format ( '  New tetra: ',4i7)

  return
end
