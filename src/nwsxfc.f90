subroutine nwsxfc ( k, i, ifac, npt, sizht, bf_num, nfc, bf_max, fc_max, &
  bf, fc, ht, nsmplx, hdavbf, hdavfc, front, back, ind, ierr )

!*****************************************************************************80
!
!! NWSXFC creates new simplices from the insertion of a face vertex.
!
!  Discussion:
!
!    This routine creates new simplices in a K-D triangulation from the
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
!    Input, integer ( kind = 4 ) K, the dimension of triangulation.
!
!    Input, integer ( kind = 4 ) I, (local) index of next vertex inserted in triangulation;
!    it is assumed I is largest index so far.
!
!    Input/output, integer ( kind = 4 ) IFAC; on input, face containing vertex I is
!    FC(*,IFAC).  On output, new face with vertex I as a vertex.
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
!    Output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of queue
!    of interior faces AB...C such that AB...CI is a new simplex.
!
!    Workspace, integer IND(1:K), the local vertex indices of K-D vertices.
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
  logical              bface
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(k+4,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ifacin
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) top

  ierr = 0
  km1 = k - 1
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  front = 0
  back = 0
  top = 0
  nv = kp2
  e = fc(kp2,ifac)
  bface = ( e <= 0 )

  if ( bface ) then
    nv = kp1
  end if

  do iv = kp1, nv

    nsmplx = nsmplx + km1
    d = fc(iv,ifac)

    do j = 1, k

      ind(1:k) = fc(1:k,ifac)
      a = ind(j)
      ind(j) = d
      pos = htsrck ( k, ind, npt, sizht, fc, ht )

      if ( pos <= 0 ) then
        ierr = 400
        return
      end if

      if ( fc(kp1,pos) == a ) then
        fc(kp1,pos) = i
      else
        fc(kp2,pos) = i
      end if

      if ( 0 < fc(kp2,pos) ) then
        if ( front == 0 ) then
           front = pos
        else
          fc(kp4,back) = pos
        end if
        back = pos
      end if

      if ( msglvl == 4 ) then
        write ( *,600) fc(1:k,pos),i
      end if

      if ( iv /= kp1 ) then
        cycle
      end if

      call availk ( k, hdavfc, nfc, fc_max, fc, pos, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      ind(1:k) = fc(1:k,ifac)
      ind(j) = i
      call htinsk ( k, pos, ind, d, e, npt, sizht, fc, ht )

      if ( bface ) then
        fc(kp4,pos) = top
        top = pos
      end if

    end do

    do j = 1, km1

      do jj = j+1, k

        call availk(k,hdavfc,nfc,fc_max,fc,pos,ierr)

        if ( ierr /= 0 ) then
          return
        end if

        l = 0
        do ii = 1, j-1
          l = l + 1
          ind(l) = fc(ii,ifac)
        end do

        do ii = j+1, jj-1
          l = l + 1
          ind(l) = fc(ii,ifac)
        end do

        do ii = jj+1, k
          l = l + 1
          ind(l) = fc(ii,ifac)
        end do

        ind(km1) = d
        ind(k) = i
        a = fc(j,ifac)
        b = fc(jj,ifac)
        call htinsk(k,pos,ind,a,b,npt,sizht,fc,ht)

      end do

    end do

  end do

  if ( front /= 0 ) then
    fc(kp4,back) = 0
  end if

  call htdelk(k,ifac,npt,sizht,fc,ht)
  fc(1,ifac) = -hdavfc
  hdavfc = ifac
  ifac = pos

  if ( .not. bface ) then
    return
  end if

  ifacin = hdavfc
  e = -e
  pos = top

  do j = k, 1, -1

    if ( hdavbf /= 0 ) then
      l = hdavbf
      hdavbf = -bf(1,hdavbf)
    else
      if ( bf_max <= bf_num ) then
        ierr = 23
        return
      end if
      bf_num = bf_num + 1
      l = bf_num
    end if

    ind(j) = top
    fc(kp2,top) = -l
    nbr = bf(j,e)
    bf(k,l) = nbr
    m = -fc(kp2,nbr)

    do ii = 1, k
      if ( bf(ii,m) == ifacin ) then
        bf(ii,m) = top
        exit
      end if
    end do

    top = fc(kp4,top)

  end do

  bf(1,e) = -hdavbf
  hdavbf = e
  top = pos

  do j = k, 1, -1

    l = -fc(kp2,top)
    ii = 0

    do jj = 1, j-1
      ii = ii + 1
      bf(ii,l) = ind(jj)
    end do

    do jj = j+1, k
      ii = ii + 1
      bf(ii,l) = ind(jj)
    end do

    pos = top
    top = fc(kp4,top)
    fc(kp4,pos) = -1

  end do

  600 format ('  New simplex: ',9i7)

  return
end
