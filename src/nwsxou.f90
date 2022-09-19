subroutine nwsxou ( k, i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, &
  nsmplx, hdavbf, hdavfc, front, back, bfi, ind, ierr )

!*****************************************************************************80
!
!! NWSXOU creates new simplices for vertices outside the current convex hull.
!
!  Discussion:
!
!    This routine creates new simplices in a K-D triangulation outside
!    the convex hull by joining vertex I to visible boundary faces.
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
!    Input, integer ( kind = 4 ) FRONT, the index of front of queue (or top of stack)
!    of visible boundary faces.
!
!    Output, integer ( kind = 4 ) BACK, the index of back of queue (or bottom of stack)
!    of visible boundary faces (which become interior faces).
!
!    Output, integer ( kind = 4 ) BFI, the index of FC of a boundary face containing
!    vertex I.
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
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfnew
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
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) j
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
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) ptr
!
!  For AB...C in queue, form simplex AB...CI + K faces involving I.
!  PTR, NBR, POS are indices of FC; L, M, BFNEW indices of BF.
!
  ierr = 0
  km1 = k - 1
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  bfi = 0
  ptr = front

  do

    back = ptr
    l = -fc(kp2,ptr)
    fc(kp2,ptr) = i
    nsmplx = nsmplx + 1

    if ( msglvl == 4 ) then
      write ( *,600) (fc(j,ptr),j=1,k),i
    end if

    do e = 1, k

      ind(1:k) = fc(1:k,ptr)
      d = ind(e)
      ind(e) = i
      nbr = bf(e,l)

      if ( fc(kp4,nbr) /= -1 ) then
        if ( fc(kp2,nbr) == i ) then
          cycle
        end if
      end if

      call availk(k,hdavfc,nfc,fc_max,fc,pos,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      m = -fc(kp2,nbr)

      do ii = 1, k
        if ( bf(ii,m) == ptr ) then
          j = ii
          exit
        end if
      end do

      if ( fc(kp4,nbr) /= -1 ) then

        call htinsk(k,pos,ind,d,fc(j,nbr),npt,sizht,fc,ht)

      else

        if ( hdavbf /= 0 ) then

          bfnew = hdavbf
          hdavbf = -bf(1,hdavbf)

        else

          if ( bf_max <= bf_num ) then
            ierr = 23
            return
          end if

          bf_num = bf_num + 1
          bfnew = bf_num

        end if

        call htinsk(k,pos,ind,d,-bfnew,npt,sizht,fc,ht)
        fc(kp4,pos) = bfi
        bfi = pos
        bf(j,m) = pos
        bf(k,bfnew) = nbr
        bf(1:km1,bfnew) = 0

      end if

    end do

    bf(1,l) = -hdavbf
    hdavbf = l
    ptr = fc(kp4,ptr)

    if ( ptr == 0 ) then
      exit
    end if

  end do
!
!  Set BF(1:K-1,*) fields for new boundary faces, which are in stack
!  with top pointer BFI.
!
  ptr = bfi

70 continue

  l = -fc(kp2,ptr)

  do e = 1, km1

    if ( 0 < bf(e,l) ) then
      cycle
    end if

    a = fc(e,ptr)
    d = fc(kp1,ptr)

80  continue

    ind(1:k) = fc(1:k,ptr)
    ind(e) = d
    nbr = htsrck(k,ind,npt,sizht,fc,ht)

    if ( nbr <= 0 ) then
      ierr = 400
      return
    end if

    if ( 0 < fc(kp2,nbr) ) then

      if ( fc(kp1,nbr) == a ) then
        a = d
        d = fc(kp2,nbr)
      else
        a = d
        d = fc(kp1,nbr)
      end if

      go to 80

    end if

    bf(e,l) = nbr
    m = -fc(kp2,nbr)

    do ii = 1, km1
      if ( fc(ii,nbr) == d ) then
        bf(ii,m) = ptr
        exit
      end if
    end do

  end do

  pos = ptr
  ptr = fc(kp4,ptr)
  fc(kp4,pos) = -1

  if ( ptr /= 0 ) then
    go to 70
  end if

  600 format ( '  New simplex: ',9i7)

  return
end
