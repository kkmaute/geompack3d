subroutine swaphs ( k, i, npt, sizht, bf_num, nfc, bf_max, fc_max, vcl, vm, bf, &
  fc, ht, nsmplx, hdavbf, hdavfc, front, back, ifac, bfi, ind, indf, mv, loc, &
  zpn, alpha, mat, ierr )

!*****************************************************************************80
!
!! SWAPHS swaps faces in a KD triangulation.
!
!  Discussion:
!
!    This routine swaps faces (apply local transformations) in a
!    K-D triangulation based on the empty hypercircumsphere criterion
!    until all faces are locally optimal, where I is index of new vertex
!    added to triangulation.
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
!    Input, integer ( kind = 4 ) I, the local index of next vertex inserted in triangulation;
!    it is assumed I is largest index so far.
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
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of queue
!    of interior faces for which hypersphere test is applied.
!
!    Output, integer ( kind = 4 ) IFAC, the index of last face for which sphere test
!    applied (if swap then index of new interior face), or 0.
!
!    Output, integer ( kind = 4 ) BFI, the index of FC of a boundary face containing vertex
!    I if a degenerate local transformation is applied, or 0.
!
!    Workspace, integer IND(1:K+1), the indices in VCL or local vertex indices.
!
!    Workspace, integer INDF(1:K+1), the pivot indices or face indices.
!
!    Workspace, integer MV(1:K+1), the missing local vertex indices used in
!    degenerate cases.
!
!    Workspace, integer LOC(1:K), the used by routine SWAPDG.
!
!    Workspace, integer ZPN(1:K), the used by routine SWAPDG.
!
!    Workspace, real ( kind = 8 ) ALPHA(1:K+1), the barycentric coordinates
!    or hypersphere center.
!
!    Workspace, real ( kind = 8 ) MAT(1:K,1:K), the matrix used for solving
!    system of linear equations.
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
  real    ( kind = 8 ) alpha(k+1)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(k,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfn
  integer ( kind = 4 ) bfp
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  logical              degen
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
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind(k+1)
  integer ( kind = 4 ) indf(k+1)
  integer ( kind = 4 ) izero
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jn
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) loc(k)
  real    ( kind = 8 ) mat(k,k)
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) mv(k+1)
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) num
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) ptrn
  real    ( kind = 8 ) radsq
  integer ( kind = 4 ) sneg
  integer ( kind = 4 ) spos
  integer ( kind = 4 ) szero
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(k,*)
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vm(npt)
  integer ( kind = 4 ) zpn(k)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  km1 = k - 1
  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  ifac = 0
  bfi = 0

10 continue

  if ( front == 0 ) then
    return
  end if

  pos = front
  front = fc(kp4,pos)

  if ( fc(2,pos) == 0 ) then
    fc(1,pos) = -hdavfc
    hdavfc = pos
    go to 10
  end if

  ifac = pos
  fc(kp4,pos) = -1
  ind(1:k) = vm(fc(1:k,pos))
  d = fc(kp1,pos)
  e = fc(kp2,pos)

  if ( d == i ) then
    d = e
    e = i
  end if

  ind(kp1) = vm(d)
  ve = vm(e)

  if ( msglvl == 4 ) then
    write ( *,600) pos,(fc(ii,pos),ii=1,k),d,e
  end if

  call ccsphk(k,.true.,ind,vcl,vcl(1,ve),alpha,radsq,in,mat,indf)

  if ( in < 1 ) go to 10

  call baryck(k,ind,vcl,vcl(1,ve),alpha,degen,mat,indf)

  if ( degen ) then
    ierr = 401
    return
  end if

  sneg = 2
  szero = 0

  do j = 1, k

    if ( alpha(j) < -tol ) then
      alpha(j) = -1.0D+00
      sneg = sneg + 1
    else if ( alpha(j) <= tol ) then
      alpha(j) = 0.0D+00
      szero = szero + 1
      izero = j
    else
      alpha(j) = 1.0D+00
    end if

  end do

  spos = kp2 - sneg - szero

  if ( msglvl == 4 ) then
    write ( *,610) (alpha(j),j=1,k)
  end if

  if ( sneg < 2 .or. spos < 2 ) then

    ierr = 405
    return
!
!  Swap SNEG simplices for SPOS if possible.  Process old interior,
!  then new interior, then boundary faces. E = I.
!
  else if ( szero == 0 ) then

    if ( 3 <= sneg ) then

      jj = 0

      do j = 1, k

        if ( alpha(j) /= -1.0D+00 ) then
          cycle
        end if

        ind(1:k) = fc(1:k,pos)
        ind(j) = d
        ptr = htsrck(k,ind,npt,sizht,fc,ht)
        if ( ptr <= 0) go to 450
        jj = jj + 1
        indf(jj) = ptr

        if ( fc(kp1,ptr) /= e .and. fc(kp2,ptr) /= e ) then
          if ( msglvl == 4 ) then
            write ( *,620) sneg,spos
          end if
          go to 10
        end if

      end do

      jj = 0

      do j = 1, k

        if ( alpha(j) /= -1.0D+00 ) then
          cycle
        end if

        jj = jj + 1
        ptr = indf(jj)
        call htdelk(k,ptr,npt,sizht,fc,ht)

        if ( 0 <= fc(kp4,ptr) ) then
          fc(2,ptr) = 0
        else
          fc(1,ptr) = -hdavfc
          hdavfc = ptr
        end if

        ind(1:k) = fc(1:k,pos)
        ind(j) = e
        call htsdlk(k,ind,npt,sizht,fc,ht,ptr)
        if ( ptr <= 0) go to 450
        fc(1,ptr) = -hdavfc
        hdavfc = ptr

      end do

    end if

    num = 1
    nsmplx = nsmplx + (spos - sneg)

    if ( msglvl == 4 ) then
      write ( *,630) sneg,spos
    end if
!
!  Swap SNEG simplices for SPOS if possible. Two simultaneous
!  local transformations may be needed. Process old and new
!  faces of degenerate local transformation, then other old
!  interior, other new interior, other boundary faces. E = I.
!
  else if ( szero == 1 ) then

    if ( 3 <= sneg ) then

      jj = 0

      do j = 1,k

        if ( alpha(j) /= -1.0D+00 ) then
          cycle
        end if

        ind(1:k) = fc(1:k,pos)
        ind(izero) = d
        ind(j) = e
        ptr = htsrck(k,ind,npt,sizht,fc,ht)

        if ( ptr <= 0 ) then
          if ( msglvl == 4 ) then
            write ( *,620) sneg,spos, 'missing face'
          end if
          go to 10
        end if

        jj = jj + 1
        indf(jj) = ptr
        mv(jj) = fc(j,pos)

      end do

    end if

    ind(1:k) = fc(1:k,pos)
    ind(izero) = e
    ptr = htsrck(k,ind,npt,sizht,fc,ht)
    if ( ptr <= 0) go to 450
    indf(sneg-1) = ptr
    ind(1:k) = fc(1:k,pos)
    a = ind(izero)
    ind(izero) = d
    ptr = htsrck(k,ind,npt,sizht,fc,ht)
    if ( ptr <= 0) go to 450
    indf(sneg) = ptr
    f = fc(kp2,ptr)

    if ( 0 < f ) then

      if ( f == a ) then
        f = fc(kp1,ptr)
      end if

      do j = 1,sneg-1

        ptr = indf(j)

        do ii = kp1,kp2

          b = fc(ii,ptr)

          if ( b /= a .and. b /= f ) then
            if ( msglvl == 4 ) then
              write ( *,620) sneg,spos, 'other vertices not common'
            end if
            go to 10
          end if

        end do

      end do

      num = 2

      do j = 1,sneg
        ptr = indf(j)
        call htdelk(k,ptr,npt,sizht,fc,ht)
        fc(1,ptr) = -hdavfc
        hdavfc = ptr
      end do

      do j = 1,k

        if ( alpha(j) /= 1.0D+00 ) then
          cycle
        end if

        call availk(k,hdavfc,nfc,fc_max,fc,ptr,ierr)

        if ( ierr /= 0 ) then
          return
        end if

        ind(1:k) = fc(1:k,pos)
        ind(izero) = d
        ind(j) = e
        call htinsk(k,ptr,ind,a,f,npt,sizht,fc,ht)

      end do

      ind(1:k) = fc(1:k,pos)
      ind(izero) = f
      call htsdlk(k,ind,npt,sizht,fc,ht,ptr)

      if ( ptr <= 0) go to 450

      if ( 0 <= fc(kp4,ptr) ) then
        fc(2,ptr) = 0
      else
        fc(1,ptr) = -hdavfc
        hdavfc = ptr
      end if

      nsmplx = nsmplx + 2*(spos - sneg)

      if ( msglvl == 4 ) then
        write ( *,630) sneg,spos, 'other vertex =',f
      end if

    else

      do j = 1,sneg-2
        b = fc(kp1,indf(j))
        if ( b /= a ) then
          if ( msglvl == 4 ) then
            write ( *,620) sneg,spos, 'other vertex not common'
          end if
          go to 10
        end if
      end do

      num = 1
      jj = sneg

      do j = 1,k

        if ( alpha(j) /= 1.0D+00 ) then
          cycle
        end if

        call availk(k,hdavfc,nfc,fc_max,fc,ptr,ierr)

        if ( ierr /= 0 ) then
          return
        end if

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

        ind(1:k) = fc(1:k,pos)
        ind(izero) = d
        ind(j) = e
        call htinsk(k,ptr,ind,a,-bfp,npt,sizht,fc,ht)
        jj = jj + 1
        indf(jj) = ptr
        mv(jj) = fc(j,pos)

      end do

      bfi = ptr
      mv(sneg-1) = d
      mv(sneg) = e

      if ( 3 <= sneg ) then

        j = sneg - 1
        ptr = indf(j)

210     continue

        if ( d < mv(j-1) ) then
          mv(j) = mv(j-1)
          indf(j) = indf(j-1)
           j = j - 1
          if ( 1 < j ) go to 210
        end if

        mv(j) = d
        indf(j) = ptr

      end if

      do j = sneg+1, kp1

        b = mv(j)
        ptr = indf(j)
        bfp = -fc(kp2,ptr)
        jn = 1
        jp = sneg + 1

        if ( jp == j ) then
          jp = jp + 1
        end if

        do jj = 1, k

          c = fc(jj,ptr)

          if ( c == mv(jp) ) then

            bf(jj,bfp) = indf(jp)
            jp = jp + 1

            if ( jp == j ) then
              jp = jp + 1
            end if

            if ( kp1 < jp ) then
              jp = kp1
            end if

          else

            ptrn = indf(jn)
            bfn = -fc(kp2,ptrn)
            ii = 1

            do while (fc(ii,ptrn) /= b)
              ii = ii + 1
            end do

            nbr = bf(ii,bfn)
            bf(jj,bfp) = nbr
            bfn = -fc(kp2,nbr)
            ii = 1

            do while (bf(ii,bfn) /= ptrn)
              ii = ii + 1
            end do

            bf(ii,bfn) = ptr
            jn = jn + 1

          end if

        end do

      end do

      do j = 1,sneg
        ptr = indf(j)
        bfp = -fc(kp2,ptr)
        call htdelk(k,ptr,npt,sizht,fc,ht)
        fc(1,ptr) = -hdavfc
        hdavfc = ptr
        bf(1,bfp) = -hdavbf
        hdavbf = bfp
      end do

      nsmplx = nsmplx + (spos - sneg)

      if ( msglvl == 4 ) then
        write ( *,630) sneg,spos, 'boundary degenerate'
      end if

    end if

    if ( 3 <= sneg ) then

      do j = 1, k

        if ( alpha(j) /= -1.0D+00 ) then
          cycle
        end if

        do in = 1, num+num

          ind(1:k) = fc(1:k,pos)

          if ( in == 1 .or. in == 3 ) then
            ind(j) = d
          else
            ind(j) = e
          end if

          if ( 3 <= in ) then
            ind(izero) = f
          end if

          call htsdlk(k,ind,npt,sizht,fc,ht,ptr)
          if ( ptr <= 0) go to 450

          if ( 0 <= fc(kp4,ptr) ) then
            fc(2,ptr) = 0
          else
            fc(1,ptr) = -hdavfc
            hdavfc = ptr
          end if

        end do

      end do

    end if
!
!  Call SWAPDG for case of 2 <= SZERO.
!
  else

    call swapdg(k,pos,d,i,sneg,spos,szero,alpha,npt,sizht,bf_num, &
      nfc,bf_max,fc_max,vcl,vm,bf,fc,ht,nsmplx,hdavbf,hdavfc, &
      front,back,ifac,bfi,ind,indf,mv,loc,zpn, ierr )

      if ( ierr /= 0 ) then
        return
      end if

    go to 10

  end if
!
!  Common code for cases of SZERO = 0 and 1.
!
  if ( 4 <= sneg ) then

    do j = 1,km1

      if ( alpha(j) /= -1.0D+00 ) then
        cycle
      end if

      do jj = j+1,k

        if ( alpha(jj) /= -1.0D+00 ) then
          cycle
        end if

        do in = 1,num

          ind(1:k) = fc(1:k,pos)
          ind(j) = d
          ind(jj) = e
          if ( in == 2) ind(izero) = f
          call htsdlk(k,ind,npt,sizht,fc,ht,ptr)
          if ( ptr <= 0) go to 450
          fc(1,ptr) = -hdavfc
          hdavfc = ptr

        end do

      end do

    end do

  end if

  do j = 1,k

    if ( alpha(j) /= 1.0D+00 ) then
      cycle
    end if

    a = fc(j,pos)

    do jj = j+1, k

      if ( alpha(jj) /= 1.0D+00 ) then
        cycle
      end if

      b = fc(jj,pos)

      do in = 1,num

        call availk ( k, hdavfc, nfc, fc_max, fc, ptr, ierr )

        if ( ierr /= 0 ) then
          return
        end if

        ind(1:k) = fc(1:k,pos)
        ind(j) = d
        ind(jj) = e
        if ( in == 2) ind(izero) = f
        call htinsk(k,ptr,ind,a,b,npt,sizht,fc,ht)

      end do

    end do

  end do

  ifac = ptr

  do j = 1, k

    if ( alpha(j) /= 1.0D+00 ) then
      cycle
    end if

    a = fc(j,pos)

    do in = 1, num+num

      ind(1:k) = fc(1:k,pos)

      if ( in == 1 .or. in == 3 ) then
        ind(j) = d
        b = e
      else
        ind(j) = e
        b = d
      end if

      if ( 3 <= in ) then
        ind(izero) = f
      end if

      call updatk(k,ind,a,b,i,npt,sizht,front,back,fc,ht, ierr )

      if ( ierr /= 0 ) then
        return
      end if

    end do

  end do

  if ( 3 <= sneg ) then

    do j = 1,k

      if ( alpha(j) /= 1.0D+00 ) then
        cycle
      end if

      a = fc(j,pos)

      do jj = 1, k

        if ( alpha(jj) /= -1.0D+00 ) then
          cycle
        end if

        b = fc(jj,pos)

        do in = 1, num
          ind(1:k) = fc(1:k,pos)
          ind(j) = d
          ind(jj) = e
          if ( in == 2) ind(izero) = f
          call updatk(k,ind,a,b,i,npt,sizht,front,back,fc,ht, ierr )

          if ( ierr /= 0 ) then
            return
          end if

        end do

      end do

    end do

  end if

  call htdelk ( k, pos, npt, sizht, fc, ht )
  fc(1,pos) = -hdavfc
  hdavfc = pos

  go to 10

450 continue

  ierr = 400

  600 format (1x,'index=',i7,':',9i7)
  610 format (4x,'sign(1:k)= ',7f4.0)
  620 format ( '  No swap ',i2,' -',i2,3x,a)
  630 format ( '  Swap ',i2,' -',i2,3x,a,i7)

  return
end
