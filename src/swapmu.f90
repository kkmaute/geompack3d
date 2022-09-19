subroutine swapmu ( bndcon, crit, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, &
  ht, ntetra, hdavfc, front, back, ifac, ierr )

!*****************************************************************************80
!
!! SWAPMU swaps faces in a KD triangulation.
!
!  Discussion:
!
!    This routine swaps faces (apply local transformations) in a 3D
!    triangulation based on a local max-min solid angle criterion until all
!    faces are locally optimal.
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
!    Input, logical BNDCON, TRUE iff boundary faces are constrained (i.e. not
!    swapped by local transformations).
!
!    Input, integer ( kind = 4 ) CRIT, criterion code; 1 for (local max-min) solid angle
!    criterion, 2 for radius ratio criterion, 3 for mean ratio
!    criterion, 0 (or anything else) for no swaps.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the array of boundary face records;
!    see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of queue
!    of interior faces for which local optimality test is applied.
!
!    Output, integer ( kind = 4 ) IFAC, the index of last face processed from queue, or 0.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  real    ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) bf(3,*)
  logical              bndcon
  integer ( kind = 4 ) bfx(2)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) crit
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  logical              degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  real    ( kind = 8 ) gama
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) indx(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  real    ( kind = 8 ) m1
  real    ( kind = 8 ) m2
  real    ( kind = 8 ) m3
  real    ( kind = 8 ) m4
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nbr(2,2)
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  real    ( kind = 8 ) s(4)
  real    ( kind = 8 ) tetmu
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vf
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  ifac = 0

10 continue

  if ( front == 0 ) then
    return
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
    go to 10
  end if

  ifac = ind
  fc(7,ind) = -1
  a = fc(1,ind)
  b = fc(2,ind)
  c = fc(3,ind)
  d = fc(4,ind)
  e = fc(5,ind)
  va = vm(a)
  vb = vm(b)
  vc = vm(c)
  vd = vm(d)
  ve = vm(e)

  if ( msglvl == 4 ) then
    write ( *,600) ind,a,b,c,d,e
  end if

  call baryth(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve), &
    alpha,degen)

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

  if ( kneg == 1 .and. kzero == 0 ) then

    m1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
    m2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
    beta = min ( m1, m2 )
    m1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vd),vcl(1,ve),s)
    m2 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
    m3 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
    gama = min ( m1, m2, m3 )
!
!  Swap 2 tetrahedra for 3.
!
    if ( beta + tol < gama ) then

      call updatf(a,b,d,c,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,b,e,c,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
      call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htdel(ind,npt,sizht,fc,ht)
      call htins(ind,a,d,e,b,c,npt,sizht,fc,ht)
      call availf ( hdavfc, nfc, fc_max, fc, ind, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      call htins(ind,b,d,e,a,c,npt,sizht,fc,ht)
      call availf(hdavfc,nfc,fc_max,fc,ind,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
      ntetra = ntetra + 1

      if ( msglvl == 4 ) then
        write ( *,610)
      end if

    end if
!
!  Relabel so edge AB would be deleted by swap.
!
  else if ( kneg == 2 .and. kzero == 0 ) then

    if ( alpha(1) < 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) < 0.0D+00 ) then
      call i4_swap ( b, c )
    end if

    ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
    if ( ind1 <= 0 ) go to 50

    if ( fc(4,ind1) == e .or. fc(5,ind1) == e ) then

      va = vm(a)
      vb = vm(b)
      vc = vm(c)
      m1 =tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
      m2 =tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
      m3 =tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vd),vcl(1,ve),s)
      beta = min ( m1, m2, m3 )
      m1 =tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
      m2 =tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
      gama = min ( m1, m2 )
!
!  Swap 3 tetrahedra for 2.
!
      if ( beta + tol < gama ) then

        call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(a,d,e,b,c,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,d,e,a,c,0,npt,sizht,front,back,fc,ht, ierr )

        if ( ierr /= 0 ) then
          return
        end if

        call htdel(ind,npt,sizht,fc,ht)
        call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
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

        ind1 = htsrc(a,b,e,npt,sizht,fc,ht)

        if ( ind1 <= 0 ) then
          go to 50
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

        if ( msglvl == 4 ) then
          write ( *,620) c,d,e
        end if

      end if

    end if
!
!  Relabel vertices so that DE intersects AB.
!  Also swap if necessary to make A < B and D < E.
!
  else if ( kneg == 1 .and. kzero == 1 ) then

    if ( alpha(1) == 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) == 0.0D+00 ) then
      call i4_swap ( b, c )
    end if

    if ( b < a ) then
      call i4_swap ( a, b )
    end if

    if ( e < d ) then
      call i4_swap ( d, e )
    end if

    ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
    ind2 = htsrc(a,b,e,npt,sizht,fc,ht)

    if ( ind1 <= 0 .or. ind2 <= 0 ) then
      go to 50
    end if

    if ( fc(4,ind1) == c ) then
      f = fc(5,ind1)
    else
      f = fc(4,ind1)
    end if

    if ( fc(4,ind2) == c ) then
      g = fc(5,ind2)
    else
      g = fc(4,ind2)
    end if

    if ( f <= 0 .and. g <= 0 ) then

      if ( .not. bndcon ) then

        va = vm(a)
        vb = vm(b)
        vc = vm(c)
        vd = vm(d)
        ve = vm(e)

        m1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
        m2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
        beta = min ( m1, m2 )
        m1 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
        m2 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
        gama = min ( m1, m2 )
!
!  Swap 2 tetrahedra for 2.
!
        if ( beta + tol < gama ) then

          call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )

          if ( ierr /= 0 ) then
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)
          call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
          call htdel(ind1,npt,sizht,fc,ht)
          call htins(ind1,a,d,e,c,fc(5,ind1),npt,sizht,fc,ht)
          call htdel(ind2,npt,sizht,fc,ht)
          call htins(ind2,b,d,e,c,fc(5,ind2),npt,sizht,fc,ht)
          indx(1) = ind1
          indx(2) = ind2
          bfx(1) = -fc(5,ind1)
          bfx(2) = -fc(5,ind2)
          dd = d

          do j = 1, 2

            if ( j == 2 ) then
              dd = e
            end if

            if ( dd < a ) then
              nbr(j,1) = bf(3,bfx(j))
              nbr(j,2) = bf(2,bfx(j))
            else if ( dd < b ) then
              nbr(j,1) = bf(3,bfx(j))
              nbr(j,2) = bf(1,bfx(j))
            else
              nbr(j,1) = bf(2,bfx(j))
              nbr(j,2) = bf(1,bfx(j))
            end if

          end do

          aa = a
          k = -fc(5,nbr(1,2))

          do j = 1, 2

            if ( j == 2 ) then
              aa = b
              k = -fc(5,nbr(2,1))
            end if

            if ( aa < d ) then
              bf(1,bfx(j)) = indx(3-j)
              bf(2,bfx(j)) = nbr(2,j)
              bf(3,bfx(j)) = nbr(1,j)
            else if ( aa < e ) then
              bf(1,bfx(j)) = nbr(2,j)
              bf(2,bfx(j)) = indx(3-j)
              bf(3,bfx(j)) = nbr(1,j)
            else
              bf(1,bfx(j)) = nbr(2,j)
              bf(2,bfx(j)) = nbr(1,j)
              bf(3,bfx(j)) = indx(3-j)
            end if

            if ( bf(1,k) == indx(j) ) then
              bf(1,k) = indx(3-j)
            else if ( bf(2,k) == indx(j) ) then
              bf(2,k) = indx(3-j)
            else
              bf(3,k) = indx(3-j)
            end if

          end do

          if ( msglvl == 4 ) then
            write ( *,630) a,b,d,e
          end if

        end if

      end if

    else if ( f == g ) then

      va = vm(a)
      vb = vm(b)
      vc = vm(c)
      vd = vm(d)
      ve = vm(e)
      vf = vm(f)
      m1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
      m2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
      m3 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vd),vcl(1,vf),s)
      m4 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,ve),vcl(1,vf),s)
      beta = min ( m1, m2, m3, m4 )
      m1 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
      m2 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
      m3 = tetmu(crit,vcl(1,va),vcl(1,vd),vcl(1,ve),vcl(1,vf),s)
      m4 = tetmu(crit,vcl(1,vb),vcl(1,vd),vcl(1,ve),vcl(1,vf),s)
      gama = min ( m1, m2, m3, m4 )
!
!  Swap 4 tetrahedra for 4.
!
      if ( beta + tol < gama ) then

        call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(a,d,f,b,e,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(a,e,f,b,d,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,d,f,a,e,0,npt,sizht,front,back,fc,ht, ierr )
        call updatf(b,e,f,a,d,0,npt,sizht,front,back,fc,ht, ierr )

        if ( ierr /= 0 ) then
          return
        end if

        call htdel(ind,npt,sizht,fc,ht)
        call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
        ind = htsrc(a,b,f,npt,sizht,fc,ht)
        if ( ind <= 0) go to 50
        call htdel(ind,npt,sizht,fc,ht)

        if ( 0 <= fc(7,ind) ) then
          fc(2,ind) = 0
          call availf(hdavfc,nfc,fc_max,fc,ind,ierr)

          if ( ierr /= 0 ) then
            return
          end if

        end if

        call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)
        call htdel(ind1,npt,sizht,fc,ht)
        j = fc(7,ind1)
        call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)
        fc(7,ind1) = j
        call htdel(ind2,npt,sizht,fc,ht)
        j = fc(7,ind2)
        call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)
        fc(7,ind2) = j

        if ( fc(7,ind1) == -1 ) then

          fc(7,ind1) = 0

          if ( front == 0 ) then
            front = ind1
          else
            fc(7,back) = ind1
          end if

          back = ind1

        end if

        if ( fc(7,ind2) == -1 ) then

          fc(7,ind2) = 0

          if ( front == 0 ) then
            front = ind2
          else
            fc(7,back) = ind2
          end if

          back = ind2

        end if

        if ( msglvl == 4 ) then
          write ( *,640) a,b,d,e,f
        end if

      end if

    end if

  end if

  go to 10

50 continue

  ierr = 300

  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 2-2: edge ',2i7,' repl by ',2i7)
  640 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)

  return
end
