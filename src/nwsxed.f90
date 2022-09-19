subroutine nwsxed ( k, i, ifac, nv, indf, npt, sizht, bf_num, nfc, bf_max, &
  fc_max, bf, fc, ht, nsmplx, hdavbf, hdavfc, front, back, ind, loc, ierr )

!*****************************************************************************80
!
!! NWSXED creates new simplices from insertion of an interior vertex.
!
!  Discussion:
!
!    This routine creates new simplices in a K-D triangulation from insertion
!    of vertex I in interior of (NV-1)-facet of FC(*,IFAC).
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
!    triangulation; it is assumed I is largest index so far.
!
!    Input/output, integer ( kind = 4 ) IFAC, a face containing vertex I is FC(*,IFAC).
!    On output, a new face with vertex I as a vertex.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices in facet containing I,
!    2 <= NV <= K-1.
!
!    Input, integer ( kind = 4 ) INDF(1:NV), the local indices of facet vertices
!    in increasing order.
!
!    Input, integer ( kind = 4 ) NPT, the number of K-D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, BF_NUM, the number of positions used in BF array.
!
!    Input, integer ( kind = 4 ) NFC, the number of positions used in FC array.
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
!    Workspace, integer LOC(1:K), a permutation of 1 to K.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,bf_max)
  logical              bface
  integer ( kind = 4 ) bot
  integer ( kind = 4 ) botn
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
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) indf(nv)
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) kv
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loc(k)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nbrn
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) nvm1
  integer ( kind = 4 ) nvp1
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) top
  integer ( kind = 4 ) topb
  integer ( kind = 4 ) topn
!
!  TOP, BOT, PTR are top, bottom, current pointers to list of faces
!  containing given (NV-1)-facet. TOPN, BOTN are top, bottom ptrs to
!  list of new boundary faces, in same relative order as other list.
!
  ierr = 0
  km1 = k - 1
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  nvm1 = nv - 1
  nvp1 = nv + 1
  front = 0
  back = 0
  top = ifac
  bot = top
  ptr = top
  fc(kp4,top) = 0
  topn = 0

10 continue

  j = 0
  jj = nv
  l = 1

  do ii = 1, k
    if ( fc(ii,ptr) == indf(l) ) then
      j = j + 1
      loc(j) = ii
      if ( l < nv ) then
        l = l + 1
      end if
    else
      jj = jj + 1
      loc(jj) = ii
    end if
  end do

  d = fc(kp1,ptr)
  e = fc(kp2,ptr)
  bface = ( e <= 0 )

  if ( bface ) then
    kv = kp1
  else
    kv = kp2
  end if

  do j = 1, nv

    call availk ( k, hdavfc, nfc, fc_max, fc, pos, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    ind(1:k) = fc(1:k,ptr)
    ind(loc(j)) = i
    call htinsk(k,pos,ind,d,e,npt,sizht,fc,ht)
    ifac = pos

    if ( bface ) then
      if ( topn == 0 ) then
        topn = pos
      else
        fc(kp4,botn) = pos
      end if
      botn = pos
    end if

  end do

  do iv = kp1, kv

    d = fc(iv,ptr)

    do j = 1, nv

      ind(1:k) = fc(1:k,ptr)
      a = ind(loc(j))
      ind(loc(j)) = d
      pos = htsrck ( k, ind, npt, sizht, fc, ht )

      if ( pos <= 0 ) then
        ierr = 400
        return
      end if

      if ( j == 1 ) then
        if ( fc(kp1,pos) == i .or. fc(kp2,pos) == i ) then
          go to 100
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
        write ( *,600) (fc(ii,pos),ii=1,k),i
      end if

    end do

    do j = 1, nvm1
      do jj = j+1, nv

        call availk(k,hdavfc,nfc,fc_max,fc,pos,ierr)

        if ( ierr /= 0 ) then
          return
        end if

        ind(1:k) = fc(1:k,ptr)
        a = ind(loc(j))
        ind(loc(j)) = d
        b = ind(loc(jj))
        ind(loc(jj)) = i
        call htinsk(k,pos,ind,a,b,npt,sizht,fc,ht)

      end do
    end do

    nsmplx = nsmplx + nvm1

100 continue

    do j = nvp1,k

      ind(1:k) = fc(1:k,ptr)
      ind(loc(j)) = d
      pos = htsrck(k,ind,npt,sizht,fc,ht)

      if ( pos <= 0 ) then
        ierr = 400
        return
      end if

      if ( fc(kp4,pos) == -1 ) then
        fc(kp4,bot) = pos
        fc(kp4,pos) = 0
        bot = pos
      end if

    end do

  end do

  ptr = fc(kp4,ptr)

  if ( ptr /= 0 ) then
    go to 10
  end if

  if ( front /= 0 ) then
    fc(kp4,back) = 0
  end if

  if ( topn /= 0 ) then
    fc(kp4,botn) = 0
  end if
!
!  Delete faces in TOP-BOT list and process boundary faces.
!
  topb = 0

140 continue

  ptr = top
  top = fc(kp4,top)
  e = fc(kp2,ptr)
  call htdelk ( k, ptr, npt, sizht, fc, ht )

  if ( 0 < e ) then

    fc(1,ptr) = -hdavfc
    hdavfc = ptr

  else

    e = -e
    fc(kp2,ptr) = 0
    fc(kp4,ptr) = topb
    topb = ptr
    j = 0
    jj = nv
    l = 1

    do ii = 1, k
      if ( fc(ii,ptr) == indf(l) ) then
        j = j + 1
        loc(j) = ii
        if ( l < nv ) then
          l = l + 1
        end if
      else
        jj = jj + 1
        loc(jj) = ii
      end if
    end do

    pos = topn

    do j = 1, nv

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

      ind(j) = topn
      fc(loc(j),ptr) = -topn
      fc(kp2,topn) = -l
      nbr = bf(loc(j),e)
      bf(k,l) = nbr
      m = -fc(kp2,nbr)

      do ii = 1, k
        if ( bf(ii,m) == ptr ) then
          bf(ii,m) = topn
          exit
        end if
      end do

      topn = fc(kp4,topn)

    end do

    topn = pos

    do j = 1, nv

      l = -fc(kp2,topn)
      iv = nv
      jj = 1
      if ( j == jj ) then
        jj = 2
      end if

      do ii = 1, km1

        if ( fc(ii,topn) == indf(jj) ) then

          bf(ii,l) = ind(jj)

          if ( jj < nv ) then
            jj = jj + 1
            if ( j == jj .and. jj < nv ) then
              jj = jj + 1
            end if
          end if

        else

          iv = iv + 1
          nbr = bf(loc(iv),e)

          if ( fc(kp2,nbr) < 0 ) then

            bf(ii,l) = nbr

          else

            a = 0

            do b = 1, k
              nbrn = fc(b,nbr)
              if ( nbrn < 0 ) then
                a = a + 1
                if ( a == j ) then
                  exit
                end if
              end if
            end do

            nbrn = -nbrn
            bf(ii,l) = nbrn
            m = -fc(kp2,nbrn)

            do b = 1, k
              if ( bf(b,m) == ptr ) then
                bf(b,m) = topn
                go to 220
              end if
            end do

          end if

        end if

220     continue

      end do

      pos = topn
      topn = fc(kp4,topn)
      fc(kp4,pos) = -1

    end do

    bf(1,e) = -hdavbf
    hdavbf = e

  end if

  if ( top /= 0 ) then
    go to 140
  end if

  do while ( 0 < topb )
    ptr = topb
    topb = fc(kp4,topb)
    fc(1,ptr) = -hdavfc
    hdavfc = ptr
  end do

  600 format ( '  New simplex: ',9i7)

  return
end
