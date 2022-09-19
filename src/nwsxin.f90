subroutine nwsxin ( k, i, ifac, ivrt, npt, sizht, nfc, fc_max, fc, ht, nsmplx, &
  hdavfc, front, back, ind, ierr )

!*****************************************************************************80
!
!! NWSXIN creates new simplices from the insertion of an interior vertex.
!
!  Discussion:
!
!    This routine creates new simplices in a K-D triangulation from the
!    insertion of vertex I in interior of simplex.
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
!    Input, integer ( kind = 4 ) I, the (local) index of next vertex inserted in
!    triangulation.
!
!    Input, integer ( kind = 4 ) IFAC, the face of simplex containing vertex I is FC(*,IFAC).
!
!    Input, integer ( kind = 4 ) IVRT, the K+1 or K+2 where K+1st vertex of simplex
!    = FC(IVRT,IFAC).
!
!    Input, integer ( kind = 4 ) NPT, the number of K-D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:FC_MAX), the array of face records; see
!    routine DTRISK.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NSMPLX, the number of simplices in triangulation.
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
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) back
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(k+4,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) pos

  ierr = 0
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  front = 0
  back = 0
  nsmplx = nsmplx + k
  d = fc(ivrt,ifac)

  do j = 0,k

    ind(1:k) = fc(1:k,ifac)

    if ( j == 0 ) then
      a = d
      pos = ifac
    else
      a = ind(j)
      ind(j) = d
      pos = htsrck(k,ind,npt,sizht,fc,ht)
      if ( pos <= 0 ) then
        ierr = 400
        return
      end if
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

  end do

  if ( front /= 0 ) then
    fc(kp4,back) = 0
  end if

  a = fc(kp1,ifac)
  fc(kp1,ifac) = d

  do j = 1, k

    do jj = j+1, kp1

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

      do ii = jj+1, kp1
        l = l + 1
        ind(l) = fc(ii,ifac)
      end do

      ind(k) = i
      d = fc(j,ifac)
      e = fc(jj,ifac)
      call htinsk(k,pos,ind,d,e,npt,sizht,fc,ht)

    end do

  end do

  fc(kp1,ifac) = a

  600 format ( '  New simplex: ',9i7)

  return
end
