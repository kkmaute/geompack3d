subroutine swapdg ( k, pos, d, i, sneg, spos, szero, alpha, npt, sizht, bf_num, &
  nfc, bf_max, fc_max, vcl, vm, bf, fc, ht, nsmplx, hdavbf, hdavfc, front, &
  back, ifac, bfi, ind, indf, mv, loc, zpn, ierr )

!*****************************************************************************80
!
!! SWAPDG applies swaps in a KD triangulation.
!
!  Discussion:
!
!    This routine applies simultaneous degenerate local transformations or
!    swaps in a K-D triangulation, where a swap is applied to facets
!    of dimension <= K-2.
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
!    Input, integer ( kind = 4 ) POS, the position of face in FC for which swap is
!    to be applied.
!
!    Input, integer ( kind = 4 ) D, I, the local indices of vertices forming simplices
!    with face FC(*,POS); it is assumed I is largest index so far.
!
!    Input, integer ( kind = 4 ) SNEG, SPOS, SZERO, the number of negative, positive,
!    zero vertices among FC(1:K,POS) and D,I; each is at least 2.
!
!    Input, real ( kind = 8 ) ALPHA(1:K), -1.0, 1.0, or 0.0 to indicate
!    type of vertex FC(J,POS), J = 1, ..., K.
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
!    Input, real ( kind = 8 ) VCL(1:K,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
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
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of
!    queue of interior faces for which hypersphere test is applied.
!
!    Input/output, integer ( kind = 4 ) IFAC, the index of last face for which sphere
!    test applied (if swap then index of new interior face), or 0.
!
!    Input/output, integer ( kind = 4 ) BFI, the index of FC of a boundary face
!    containing vertex I if a degenerate local transformation is applied, or 0.
!
!    Workspace, integer IND(1:K), the local vertex indices.
!
!    Workspace, integer INDF(1:K), the face indices.
!
!    Workspace, integer MV(1:K+1), the missing local vertex indices.
!
!    Workspace, integer LOC(1:K), the local vertex indices of facet,
!    indices from 1 to K.
!
!    Workspace, integer ZPN(1:K), a permutation of 1 to K.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  real    ( kind = 8 ) alpha(k)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfn
  integer ( kind = 4 ) bfp
  integer ( kind = 4 ) bot(2)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dind
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
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
  integer ( kind = 4 ) indf(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jn
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) kbf(2)
  integer ( kind = 4 ) kif(2)
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) kv
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loc(k)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) mv(k+1)
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) nstrt
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) pstrt
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) ptr1
  integer ( kind = 4 ) ptrn
  integer ( kind = 4 ) sn1
  integer ( kind = 4 ) sneg
  integer ( kind = 4 ) snp
  integer ( kind = 4 ) snpm1
  integer ( kind = 4 ) spos
  integer ( kind = 4 ) spz
  integer ( kind = 4 ) spzm1
  integer ( kind = 4 ) spzm2
  integer ( kind = 4 ) szero
  integer ( kind = 4 ) top
  integer ( kind = 4 ) topd
  integer ( kind = 4 ) topn
  real    ( kind = 8 ) vcl(k,*)
  integer ( kind = 4 ) vm(npt)
  integer ( kind = 4 ) zpn(k)

  ierr = 0
  km1 = k - 1
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  sn1 = sneg + 1
  snp = sneg + spos
  snpm1 = snp - 1
  spz = spos + szero
  spzm1 = spz - 1
  spzm2 = spz - 2
  jz = 1
  pstrt = szero + 1
  jp = pstrt
  nstrt = pstrt + spos
  jn = nstrt

  do j = 1, k

    if ( alpha(j) == 0.0D+00 ) then
      zpn(jz) = j
      jz = jz + 1
    else if ( alpha(j) == 1.0D+00 ) then
      zpn(jp) = j
      jp = jp + 1
    else
      zpn(jn) = j
      jn = jn + 1
    end if

  end do
!
!  Check whether negative degenerate facets are in triangulation,
!  and put their face indices in INDF.
!
  a = fc(zpn(1),pos)

  if ( 3 <= sneg ) then

    jj = 2

    do jn = nstrt, k

      ind(1:k) = fc(1:k,pos)
      ind(zpn(1)) = d
      ind(zpn(jn)) = i
      ptr = htsrck(k,ind,npt,sizht,fc,ht)

      if ( ptr <= 0 ) then
        if ( msglvl == 4 ) then
          write ( *,620) sneg,spos, 'dg: missing facet'
        end if
        return
      else if ( fc(kp1,ptr) /= a .and. fc(kp2,ptr) /= a ) then
        if ( msglvl == 4 ) then
          write ( *,620) sneg,spos,'dg: other vertex not common'
        end if
        return
      end if

      jj = jj + 1
      indf(jj) = ptr

    end do

  end if

  ind(1:k) = fc(1:k,pos)

  ind(zpn(1)) = i
  ptr = htsrck(k,ind,npt,sizht,fc,ht)
  if ( ptr <= 0) go to 550

  indf(1) = ptr
  ind(1:k) = fc(1:k,pos)
  ind(zpn(1)) = d
  ptr = htsrck(k,ind,npt,sizht,fc,ht)
  if ( ptr <= 0) go to 550

  indf(2) = ptr
!
!  Find all faces involving the SNEG negative facets.  None of these
!  faces are in queue.  These faces are kept in 2 lists, first for
!  those containing first facet, second for others.
!
  bot(1) = indf(1)
  bot(2) = indf(2)
  fc(kp4,bot(1)) = 0
  fc(kp4,bot(2)) = 0

  do j = 1, sneg

    l = min ( j, 2 )
    kbf(l) = 0
    kif(l) = 0
    ptr = indf(j)

    if ( 3 <= j ) then
      fc(kp4,bot(l)) = ptr
      fc(kp4,ptr) = 0
      bot(l) = ptr
    end if

    jj = 0
    jz = 2

    do ii = 1,k
      if ( fc(ii,ptr) == fc(zpn(jz),pos) ) then
        if ( jz < szero) jz = jz + 1
      else
        jj = jj + 1
        loc(jj) = fc(ii,ptr)
      end if
    end do

70  continue

    if ( 0 < fc(kp2,ptr) ) then
      kif(l) = kif(l) + 1
      kv = kp2
    else
      kbf(l) = kbf(l) + 1
      kv = kp1
    end if

    jn = 1
    jz = snpm1

    do ii = 1,k
      if ( fc(ii,ptr) == loc(jn) ) then
        if ( jn < snpm1) jn = jn + 1
      else
        jz = jz + 1
        loc(jz) = ii
      end if
    end do

    do m = kp1, kv

      e = fc(m,ptr)

      do jz = snp,k

        ind(1:k) = fc(1:k,ptr)
        ind(loc(jz)) = e
        nbr = htsrck(k,ind,npt,sizht,fc,ht)
        if ( nbr <= 0) go to 550

        if ( fc(kp4,nbr) == -1 ) then
          fc(kp4,bot(l)) = nbr
          fc(kp4,nbr) = 0
          bot(l) = nbr
        end if

      end do

    end do

    ptr = fc(kp4,ptr)
    if ( ptr /= 0) go to 70

    if ( j /= 1 ) then
      if ( kbf(2) /= kbf(1) .or. kif(2) /= kif(1) ) then
        if ( msglvl == 4 ) then
          write ( *,620) sneg,spos, 'dg: diff number of faces'
        end if
        go to 530
      end if
    end if

  end do
!
!  Check if faces containing each facet match. I = FC(K,PTR).
!
  ptr = indf(1)

130 continue

  loc(sneg) = k

  if ( 2 < sneg ) then

    jj = 1
    jn = nstrt

    do ii = 1, km1
      if ( fc(ii,ptr) == fc(zpn(jn),pos) ) then
        jj = jj + 1
        loc(jj) = ii
        jn = jn + 1
        if ( k < jn ) then
          go to 150
        end if
      end if
    end do

  end if

150 continue

  do jn = 2,sneg
    ind(1:k) = fc(1:k,ptr)
    ind(loc(jn)) = d
    nbr = htsrck(k,ind,npt,sizht,fc,ht)
    if ( nbr <= 0 ) then
      if ( msglvl == 4 ) then
        write ( *,620) sneg,spos, 'dg: unmatched face'
      end if
      go to 530
    end if
  end do

  ptr = fc(kp4,ptr)
  if ( ptr /= 0) go to 130
!
!  Apply local transformations. Process old and new faces of
!  degenerate local transformation, then other old interior,
!  other new interior, other boundary faces (for each simplex).
!
  m = (kif(1) + kif(1) + kbf(1))/szero
  nsmplx = nsmplx + m*(spos - sneg)

  if ( msglvl == 4 ) then
    write ( *,630) sneg,spos,'dg: #swaps =',m
  end if

  top = indf(2)

  do

    ptr = top
    top = fc(kp4,ptr)

    if ( 0 < fc(kp2,ptr) ) then
      call htdelk(k,ptr,npt,sizht,fc,ht)
      fc(1,ptr) = -hdavfc
      hdavfc = ptr
    end if

    if ( top == 0 ) then
      exit
    end if

  end do
!
!  TOP, TOPN are used for stacks of old, new boundary faces.
!
  top = 0
  topn = 0
  ptr = indf(1)

190 continue

  jz = 1
  jp = szero
  jn = spz
  j = pstrt

  if ( nstrt <= k ) then
    jj = nstrt
  else
    jj = pstrt
  end if

  loc(k) = k

  do ii = 1,km1

    if ( fc(ii,ptr) == fc(zpn(j),pos) ) then
      loc(jp) = ii
      jp = jp + 1
      j = j + 1
      if ( j == nstrt) j = j - 1
    else if ( fc(ii,ptr) == fc(zpn(jj),pos) ) then
      loc(jn) = ii
      jn = jn + 1
      if ( jj < k) jj = jj + 1
    else
      loc(jz) = ii
      jz = jz + 1
    end if
  end do

  e = fc(kp1,ptr)
  f = fc(kp2,ptr)

  do jp = szero,spzm1

    call availk(k,hdavfc,nfc,fc_max,fc,ptr1,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    if ( f < 0 ) then

      if ( hdavbf /= 0 ) then

        bfp = hdavbf
        hdavbf = -bf(1,hdavbf)

      else

        if ( bf_max <= bf_num ) then
          ierr = 23
          return
        else
          bf_num = bf_num + 1
          bfp = bf_num
        end if

      end if

      f = -bfp
      bfi = ptr1

    end if

    ind(1:k) = fc(1:k,ptr)
    ind(loc(jp)) = d
    call htinsk(k,ptr1,ind,e,f,npt,sizht,fc,ht)

    if ( f < 0 ) then
      fc(kp4,ptr1) = topn
      topn = ptr1
    end if

  end do

  if ( 0 < f ) then
    kv = kp2
  else
    kv = kp1
  end if

  do m = kp1,kv

    if ( m == kp2) e = f

    do jn = spz,k

      ind(1:k) = fc(1:k,ptr)
      ind(loc(jn)) = e
      call htsdlk(k,ind,npt,sizht,fc,ht,ptr1)

      if ( ptr1 <= 0 ) then
        if ( jn == spz) go to 340
        go to 550
      end if

      if ( ptr1 == pos ) then
        cycle
      end if

      if ( 0 <= fc(kp4,ptr1) ) then
        fc(2,ptr1) = 0
      else
        fc(1,ptr1) = -hdavfc
        hdavfc = ptr1
      end if

    end do

    do jn = spz,km1
      do jj = jn+1,k

        ind(1:k) = fc(1:k,ptr)
        ind(loc(jn)) = e
        ind(loc(jj)) = d
        call htsdlk(k,ind,npt,sizht,fc,ht,ptr1)
        if ( ptr1 <= 0) go to 550

        if ( ptr1 == pos ) then
          cycle
        end if

        if ( 0 <= fc(kp4,ptr1) ) then
          fc(2,ptr1) = 0
        else
          fc(1,ptr1) = -hdavfc
          hdavfc = ptr1
        end if

      end do
    end do

    do jp = szero,spzm2

      a = fc(loc(jp),ptr)

      do jj = jp+1,spzm1

        call availk(k,hdavfc,nfc,fc_max,fc,ptr1,ierr)

        if ( ierr /= 0 ) then
          return
        end if

        ind(1:k) = fc(1:k,ptr)
        ind(loc(jp)) = e
        b = ind(loc(jj))
        ind(loc(jj)) = d
        call htinsk(k,ptr1,ind,a,b,npt,sizht,fc,ht)

      end do
    end do

    ifac = ptr1

    do jp = szero,spzm1

      a = fc(loc(jp),ptr)
      b = d

      do jn = spzm1,k

        ind(1:k) = fc(1:k,ptr)
        ind(loc(jp)) = e

        if ( spzm1 < jn ) then
          b = ind(loc(jn))
          ind(loc(jn)) = d
        end if

        call updatk(k,ind,a,b,i,npt,sizht,front,back,fc,ht, ierr )

        if ( ierr /= 0 ) then
          return
        end if

      end do

    end do

340 continue

  end do

  ptr1 = ptr
  ptr = fc(kp4,ptr)

  if ( 0 < f ) then
    call htdelk(k,ptr1,npt,sizht,fc,ht)
    fc(1,ptr1) = -hdavfc
    hdavfc = ptr1
  else
    fc(kp4,ptr1) = top
    top = ptr1
  end if

  if ( ptr /= 0 ) go to 190

  if ( top == 0 ) go to 520
!
!  Set BF fields of new boundary faces.
!  Then delete old boundary faces in TOPD stack.
!
  mv(sneg) = i
  mv(sneg-1) = d

  if ( 3 <= sneg ) then

    j = 0
    do jn = nstrt,k
     j = j + 1
     mv(j) = fc(zpn(jn),pos)
    end do

    j = sneg - 1

360 continue

    if ( d < mv(j-1) ) then
      mv(j) = mv(j-1)
      j = j - 1
      if ( 1 < j ) then
        go to 360
      end if
    end if

    dind = j
    mv(j) = d
  else
    dind = 1
  end if

  j = sneg
  do jp = pstrt, nstrt-1
    j = j + 1
    mv(j) = fc(zpn(jp),pos)
  end do

  topd = 0

380 continue

  ptr = top
  top = fc(kp4,ptr)
  fc(kp4,ptr) = topd
  topd = ptr
  indf(dind) = ptr
  loc(sneg) = k

  if ( 2 < sneg ) then

    jj = 1
    jn = nstrt

    do ii = 1, km1

      if ( fc(ii,ptr) == fc(zpn(jn),pos) ) then
        jj = jj + 1
        loc(jj) = ii
        jn = jn + 1
        if ( k < jn ) go to 400
      end if

    end do

  end if

400 continue

  j = 1

  do jn = 2, sneg

    ind(1:k) = fc(1:k,ptr)
    ind(loc(jn)) = d
    nbr = htsrck(k,ind,npt,sizht,fc,ht)
    if ( nbr <= 0) go to 550

    if ( j == dind ) then
      j = j + 1
    end if

    indf(j) = nbr
    j = j + 1
    fc(kp4,nbr) = topd
    topd = nbr

  end do

  do j = snp, sn1, -1
    ptr1 = topn
    topn = fc(kp4,ptr1)
    fc(kp4,ptr1) = -1
    indf(j) = ptr1
  end do

  do j = sn1,snp

    a = mv(j)
    ptr1 = indf(j)
    bfp = -fc(kp2,ptr1)
    jn = 1

    jp = sn1
    if ( jp == j ) then
      jp = jp + 1
    end if

    jz = snp

    do jj = 1, k

      b = fc(jj,ptr1)

      if ( b == mv(jp) ) then

        bf(jj,bfp) = indf(jp)

        jp = jp + 1

        if ( jp == j ) then
          jp = jp + 1
        end if

        if ( snp < jp ) then
          jp = snp
        end if

      else if ( b == mv(jn) ) then

        ptrn = indf(jn)
        ii = 1

440     continue

        if ( fc(ii,ptrn) /= a ) then
          ii = ii + 1
          go to 440
        end if

        nbr = bf(ii,-fc(kp2,ptrn))
        bf(jj,bfp) = nbr
        bfn = -fc(kp2,nbr)
        ii = 1

450     continue

        if ( bf(ii,bfn) /= ptrn ) then
          ii = ii + 1
          go to 450
        end if

        bf(ii,bfn) = ptr1

        if ( jn < sneg ) then
          jn = jn + 1
        end if

      else

        jz = jz + 1

        if ( j == sn1 ) then

          ii = 1

460       continue

          if ( fc(ii,ptr) /= b ) then
            ii = ii + 1
            go to 460
          end if

          nbr = bf(ii,-fc(kp2,ptr))
          bfn = -fc(kp2,nbr)
          ii = 1

470       continue

          if ( bf(ii,bfn) /= ptr ) then
            ii = ii + 1
            go to 470
          end if

          f = fc(ii,nbr)
          mv(jz) = f

        else

          f = mv(jz)

        end if

        ind(1:k) = fc(1:k,ptr1)
        ind(jj) = f
        ptrn = htsrck(k,ind,npt,sizht,fc,ht)
        if ( ptrn <= 0) go to 550
        bf(jj,bfp) = ptrn

      end if

    end do

  end do

  if ( top /= 0) go to 380

510 continue

  ptr = topd
  topd = fc(kp4,ptr)
  bfp = -fc(kp2,ptr)
  call htdelk(k,ptr,npt,sizht,fc,ht)
  fc(1,ptr) = -hdavfc
  hdavfc = ptr
  bf(1,bfp) = -hdavbf
  hdavbf = bfp
  if ( topd /= 0 ) go to 510

520 continue

  fc(1,pos) = -hdavfc
  hdavfc = pos
  return

530 continue

  fc(kp4,bot(1)) = indf(2)
  top = indf(1)

540 continue

  ptr = top
  top = fc(kp4,ptr)
  fc(kp4,ptr) = -1
  if ( top /= 0) go to 540
  return

550 continue

  ierr = 400

  620 format ( '  No swap ',i2,' -',i2,3x,a)
  630 format (4x,'swap ',i2,' -',i2,3x,a,i7)

  return
end
