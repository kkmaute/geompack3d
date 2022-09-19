subroutine lfcini ( k, i, ifac, ivrt, indf, npt, sizht, bf, fc, ht, nsmplx, &
  hdavbf, hdavfc, bflag, front, back, top, ind, loc, ierr )

!*****************************************************************************80
!
!! LFCINI initializes two lists of faces.
!
!  Discussion:
!
!    This routine initializes two lists of faces and deletes some simplices,
!    faces from insertion of vertex I in interior or on boundary
!    of K-D triangulation.
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
!    Input, integer ( kind = 4 ) IFAC, the index of FC indicating simplex or face
!    containing I.
!
!    Input, integer ( kind = 4 ) IVRT, the K+1 or K+2 to indicate that FC(IVRT,IFAC)
!    is (K+1)st vertex of simplex containing I in its interior;
!    K if I lies in interior of face FC(*,IFAC); 2 to K-1 if
!    I lies in interior of facet of FC(*,IFAC) of dim IVRT-1.
!
!    Input, integer ( kind = 4 ) INDF(1:IVRT), if 2 <= IVRT <= K-1 then IVRT elements are
!    local vertex indices in increasing order of (IVRT-1)-
!    facet containing I in its interior, else not referenced.
!
!    Input, integer ( kind = 4 ) NPT, the number of K-D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, thesize of hash table HT.
!
!    Input/output, integer ( kind = 4 ) BF(1:K,1:*), the array of boundary face records;
!    see DTRISK.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:*), the array of face records;
!    see routine DTRISK.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NSMPLX, the number of simplices in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVBF, the head pointer to available BF records.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Output, logical BFLAG, TRUE iff vertex I is on boundary of triangulation.
!
!    Output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of queue of
!    interior faces that may form a new simplex with vertex I.
!
!    Output, integer ( kind = 4 ) TOP, the index of top of stack of boundary faces that
!    form a new simplex with vertex I.
!
!    Workspace, integer IND(1:K), the local vertex indices of K-D vertices.
!
!    Workspace, integer LOC(1:K), the permutation of 1 to K.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,*)
  logical              bface
  logical              bflag
  integer ( kind = 4 ) botd
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(k+4,*)
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
  integer ( kind = 4 ) indf(k-1)
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivp1
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) kbf
  integer ( kind = 4 ) kif
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) kv
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loc(k)
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) top
  integer ( kind = 4 ) topd

  ierr = 0
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  front = 0
  back = 0
  top = 0
  bflag = .false.
!
!  Vertex I is in interior of simplex.
!
  if ( kp1 <= ivrt ) then

    nsmplx = nsmplx - 1
    d = fc(ivrt,ifac)

    if ( msglvl == 4 ) then
      write ( *,600) (fc(ii,ifac),ii=1,k),d
    end if

    do j = 0, k

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

      if ( 0 <= fc(kp2,pos) ) then

        if ( front == 0 ) then
          front = pos
        else
          fc(kp4,back) = pos
        end if

        back = pos

      else

        fc(kp4,pos) = top
        top = pos
      end if

    end do
!
!  Vertex I is in interior of face.
!
  else if ( ivrt == k ) then

    e = fc(kp2,ifac)
    bflag = ( e < 0 )

    if ( bflag ) then
      kv = kp1
      nsmplx = nsmplx - 1
      e = -e
      bf(1,e) = -hdavbf
      hdavbf = e
    else
      kv = kp2
      nsmplx = nsmplx - 2
    end if

    do iv = kp1,kv

      d = fc(iv,ifac)
      if ( msglvl == 4 ) then
        write ( *,600) (fc(ii,ifac),ii=1,k),d
      end if

      do j = 1,k

        ind(1:k) = fc(1:k,ifac)
        a = ind(j)
        ind(j) = d
        pos = htsrck(k,ind,npt,sizht,fc,ht)

        if ( pos <= 0 ) then
          ierr = 400
          return
        end if

        if ( fc(kp1,pos) == a ) then
          fc(kp1,pos) = i
        else
          fc(kp2,pos) = i
        end if

        if ( 0 <= fc(kp2,pos) ) then

          if ( front == 0 ) then
            front = pos
          else
            fc(kp4,back) = pos
          end if
          back = pos
        else
          fc(kp4,pos) = top
          top = pos
        end if

      end do

    end do

    call htdelk ( k, ifac, npt, sizht, fc, ht )
    fc(1,ifac) = -hdavfc
    hdavfc = ifac
!
!  Vertex I is in interior of facet of dimension IVRT-1.
!
  else

    ivp1 = ivrt + 1
    kbf = 0
    kif = 0
    topd = ifac
    botd = topd
    ptr = topd
    fc(kp4,topd) = 0

60  continue

    j = 0
    jj = ivrt
    l = 1

    do ii = 1,k

      if ( fc(ii,ptr) == indf(l) ) then
        j = j + 1
        loc(j) = ii
        if ( l < ivrt ) then
          l = l + 1
        end if
      else
        jj = jj + 1
        loc(jj) = ii
      end if

    end do

    e = fc(kp2,ptr)
    bface = ( e < 0 )

    if ( bface ) then
      bflag = .true.
      kv = kp1
      kbf = kbf + 1
      e = -e
      bf(1,e) = -hdavbf
      hdavbf = e
    else
      kv = kp2
      kif = kif + 1
    end if

    do iv = kp1, kv

      d = fc(iv,ptr)

      do j = 1,ivrt

        ind(1:k) = fc(1:k,ptr)
        a = ind(loc(j))
        ind(loc(j)) = d
        pos = htsrck(k,ind,npt,sizht,fc,ht)

        if ( pos <= 0 ) then
          ierr = 400
          return
        end if

        if ( j == 1 ) then

          if ( fc(kp1,pos) == i .or. fc(kp2,pos) == i ) then
            go to 100
          else
            if ( msglvl == 4 ) then
              write ( *,600) (fc(ii,ptr),ii=1,k),d
            end if
          end if

        end if

        if ( fc(kp1,pos) == a ) then
          fc(kp1,pos) = i
        else
          fc(kp2,pos) = i
        end if

        if ( 0 <= fc(kp2,pos) ) then

          if ( front == 0 ) then
            front = pos
          else
            fc(kp4,back) = pos
          end if

          back = pos

        else

          fc(kp4,pos) = top
          top = pos

        end if

      end do

100   continue

      do j = ivp1,k

        ind(1:k) = fc(1:k,ptr)
        ind(loc(j)) = d
        pos = htsrck(k,ind,npt,sizht,fc,ht)

        if ( pos <= 0 ) then
          ierr = 400
          return
        end if

        if ( fc(kp4,pos) == -1 ) then
          fc(kp4,botd) = pos
          fc(kp4,pos) = 0
          botd = pos
        end if

      end do

    end do

    ptr = fc(kp4,ptr)

    if ( ptr /= 0) then
      go to 60
    end if

    nsmplx = nsmplx - ( kbf + kif + kif ) / ( kp1 - ivrt )

140 continue

    ptr = topd
    topd = fc(kp4,ptr)
    call htdelk ( k, ptr, npt, sizht, fc, ht )
    fc(1,ptr) = -hdavfc
    hdavfc = ptr

    if ( topd /= 0 ) then
      go to 140
    end if

  end if

  if ( front /= 0 ) then
    fc(kp4,back) = 0
  end if

  600 format (1x,'deleted simplex:',9i7)

  return
end
