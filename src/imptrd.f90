subroutine imptrd ( bndcon, npt, sizht, fc_max, vcl, vm, nfc, ntetra, hdavfc, &
  bf, fc, ht, ierr )

!*****************************************************************************80
!
!! IMPTRD further improves a 3D triangulation.
!
!  Discussion:
!
!    This routine further improves a given 3D triangulation towards Delaunay
!    one by using combination local transformations (not yet
!    guaranteed to produce Delaunay triangulation).
!
!    BF, FC, HT should be as output by DTRIS3 or DTRIW3,
!    except it is assumed FC(7,1:2) don't contain HDAVBF, HDAVFC.
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
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT; a good choice is a
!    prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the array of boundary face records;
!    see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) bf(3,*)
  logical              bndcon
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) cc
  real    ( kind = 8 ) center(3)
  real    ( kind = 8 ) ccradi
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  logical              first
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind0
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) ind3
  integer ( kind = 4 ) j
  real    ( kind = 8 ) maxaft
  real    ( kind = 8 ) maxbef
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mu1
  real    ( kind = 8 ) mu2
  real    ( kind = 8 ) mu3
  real    ( kind = 8 ) mu4
  real    ( kind = 8 ) mu5
  real    ( kind = 8 ) mu6
  real    ( kind = 8 ) mu7
  real    ( kind = 8 ) mu8
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  real    ( kind = 8 ) nu1
  real    ( kind = 8 ) nu2
  real    ( kind = 8 ) nu3
  real    ( kind = 8 ) nu4
  real    ( kind = 8 ) nu5
  real    ( kind = 8 ) nu6
  real    ( kind = 8 ) nu7
  logical              t1
  logical              t2
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolp1
  integer ( kind = 4 ) top
  character ( len = 3 ) type
  character ( len = 3 ) typ2
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vf
  integer ( kind = 4 ) vg
  integer ( kind = 4 ) vh
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( msglvl == 4 ) then
    write ( *,600)
  end if

  front = 0
  back = 0
  ibeg = 1
  ind = 1
  tolp1 = 1.0D+00 + tol

10 continue

  if ( fc(1,ind) <= 0 .or. fc(5,ind) <= 0 ) go to 70

  call ifacty(ind,npt,sizht,vcl,vm,fc,ht,type,a,b,c,ierr)

  if ( ierr /= 0 ) then
    return
  end if

  if ( type /= 'N32' .and. type /= 'N44' ) go to 70
  d = fc(4,ind)
  e = fc(5,ind)
  va = vm(a)
  vb = vm(b)
  vc = vm(c)
  vd = vm(d)
  ve = vm(e)

  call ccsph(.true.,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd), &
    vcl(1,ve),center,mu1,in)

  if ( in == 2 ) then
    ierr = 301
    return
  end if

  if ( in <= 0) go to 70

  if ( msglvl == 4 ) then
    write ( *,610) type,ind,a,b,c,d,e
  end if

  ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
  if ( ind1 <= 0) go to 80

  if ( fc(4,ind1) == c ) then
    f = fc(5,ind1)
  else
    f = fc(4,ind1)
  end if

  ind2 = htsrc(a,b,f,npt,sizht,fc,ht)
  if ( ind2 <= 0) go to 80
  mu1 = 0.0625D+00 / mu1
  if ( type == 'N44') go to 20
!
!  TYPE == 'N32'
!
  if ( fc(4,ind2) /= e .and. fc(5,ind2) /= e ) go to 70
  call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

  if ( msglvl == 4 ) then
    write ( *,620) aa,bb,cc,d,e,type
  end if

  if ( type /= 'T23' .and. type /= 'T32' .and. type /= 'T44' &
    .and. type /= 'N32' ) go to 70

  vf = vm(f)
  mu2 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve))
  mu3 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd))
  mu4 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,ve))
  nu1 = ccradi(vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve))
  nu2 = ccradi(vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve))
  maxbef = max ( mu1, mu2, mu3, mu4 )

  if ( type == 'T23' ) then

    nu3 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve))
    nu4 = ccradi(vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve))
    if ( max ( nu1, nu2, nu3, nu4 ) <= maxbef*tolp1) go to 70
    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  else if ( type == 'T32' ) then

    if ( cc == a ) then
      call i4_swap ( va, vb )
    end if

    mu5 = ccradi(vcl(1,vf),vcl(1,va),vcl(1,vd),vcl(1,ve))
    nu3 = ccradi(vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve))
    if ( max ( nu1, nu2, nu3 ) <= max ( maxbef, mu5 ) * tolp1 ) go to 70
    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  else if ( type == 'T44' ) then

    if ( cc == a ) then
      call i4_swap ( a, b )
      call i4_swap ( va, vb )
    end if

    ind1 = htsrc(a,d,f,npt,sizht,fc,ht)
    if ( ind1 <= 0 ) go to 80

    if ( fc(4,ind1) == b ) then
      g = fc(5,ind1)
    else
      g = fc(4,ind1)
    end if

    vg = vm(g)
    mu5 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd))
    mu6 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve))
    nu3 = ccradi(vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve))
    nu4 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve))
    nu5 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))

    if ( max ( nu1, nu2, nu3, nu4, nu5 ) &
      <= max ( maxbef, mu5, mu6 ) *tolp1 ) then
      go to 70
    end if

    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  else

    if ( cc == a ) then
      call i4_swap ( a, b )
      call i4_swap ( va, vb )
    end if

    ind1 = htsrc(a,f,d,npt,sizht,fc,ht)
    if ( ind1 <= 0 ) then
      go to 80
    end if

    if ( fc(4,ind1) == b ) then
      g = fc(5,ind1)
    else
      g = fc(4,ind1)
    end if

    ind1 = htsrc ( a, f, g, npt, sizht, fc, ht )

    if ( ind1 <= 0 ) then
      go to 80
    end if

    if ( fc(4,ind1) /= e .and. fc(5,ind1) /= e ) then
      go to 70
    end if

    call ifacty ( ind1, npt, sizht, vcl, vm, fc, ht, type, aa, bb, cc, ierr )

    if ( type /= 'T23' .and. type /= 'T32' .and. type /= 'N32' ) then
      go to 70
    end if

    vg = vm(g)
    mu5 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd))
    mu6 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve))
    nu3 = ccradi(vcl(1,vb),vcl(1,vd),vcl(1,ve),vcl(1,vf))

    if ( type == 'T23' ) then

      maxbef = max ( maxbef, mu5, mu6 )
      nu4 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve))
      nu5 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))

      if ( max ( nu1, nu2, nu3, nu4, nu5 ) <= maxbef * tolp1 ) then
        go to 70
      end if

      top = ind1

    else if ( type == 'T32' ) then

      if ( cc == a ) then
        call i4_swap ( va, vf )
      end if

      mu7 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve))
      maxbef = max ( maxbef, mu5, mu6, mu7 )
      nu4 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))

      if ( max ( nu1, nu2, nu3, nu4 ) <= maxbef * tolp1 ) then
        go to 70
      end if

      top = ind1

    else

      if ( g == c ) then
        go to 70
      end if

      if ( cc == a ) then
        call i4_swap ( a, f )
        call i4_swap ( va, vf )
      end if

      ind0 = htsrc(a,g,d,npt,sizht,fc,ht)

      if ( ind0 <= 0 ) then
        go to 80
      end if

      if ( fc(4,ind0) == f ) then
        h = fc(5,ind0)
      else
        h = fc(4,ind0)
      end if

      ind0 = htsrc(a,g,h,npt,sizht,fc,ht)

      if ( ind0 <= 0) go to 80

      if ( fc(4,ind0) /= e .and. fc(5,ind0) /= e ) go to 70

      call ifacty(ind0,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)
      if ( type /= 'T23') go to 70
      vh = vm(h)
      mu7 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,vd))
      mu8 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,ve))
      maxbef = max ( maxbef, mu5, mu6, mu7, mu8 )
      nu4 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))
      nu5 = ccradi(vcl(1,va),vcl(1,vh),vcl(1,vd),vcl(1,ve))
      nu6 = ccradi(vcl(1,vg),vcl(1,vh),vcl(1,vd),vcl(1,ve))

      if ( max ( nu1, nu2, nu3, nu4, nu5, nu6 ) <= maxbef * tolp1 ) then
        go to 70
      end if

      top = ind0
      fc(7,ind0) = ind1

    end if

    fc(7,ind1) = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  end if
!
!  TYPE == 'N44'
!
20 continue

  ind3 = ind2

  if ( fc(4,ind3) == d ) then
    g = fc(5,ind3)
  else
    g = fc(4,ind3)
  end if

  ind2 = htsrc(a,b,g,npt,sizht,fc,ht)
  if ( ind2 <= 0) go to 80

  if ( fc(4,ind2) /= e .and. fc(5,ind2) /= e ) go to 70

  call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

  if ( msglvl == 4 ) then
    write ( *,620) aa,bb,cc,e,f,type
  end if

  t1 = ( type == 'T23' .or. type == 'T32' .or. type =='T44' )

  call ifacty(ind3,npt,sizht,vcl,vm,fc,ht,typ2,aa,bb,c2,ierr)

  if ( msglvl == 4 ) then
    write ( *,620) aa,bb,c2,d,g,typ2
  end if

  t2 = (typ2 == 'T23' .or. typ2 == 'T32' .or. typ2 =='T44')
  if ( ( .not. t1 ) .and. ( .not. t2 ) ) go to 70
  first = .true.

30 continue

  if ( .not. t1 ) then

    ind2 = ind3
    type = typ2
    cc = c2
    call i4_swap ( d, e )
    call i4_swap ( vd, ve )
    call i4_swap ( f, g )

  end if

  vf = vm(f)
  vg = vm(g)

  if ( first ) then
    mu2 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve))
    mu3 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd))
    mu4 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,ve))
    mu5 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,vf))
    nu1 = ccradi(vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve))
    nu2 = ccradi(vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve))
    nu3 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve))
    nu4 = ccradi(vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve))
  else
    nu3 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve))
    nu4 = ccradi(vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve))
  end if

  if ( first ) then
    maxbef = max ( mu1, mu2, mu3, mu4, mu5 )
    maxaft = max ( nu1, nu2, nu3, nu4 )
  else
    maxaft = max ( nu1, nu2, nu3, nu4 )
  end if

  if ( type == 'T23' ) then

    nu5 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf))
    nu6 = ccradi(vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf))

    if ( max ( maxaft, nu5, nu6 ) <= maxbef * tolp1 ) then
      go to 50
    end if

    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  else if ( type == 'T32' ) then

    if ( cc == a ) then
      call i4_swap ( a, b )
      call i4_swap ( va, vb )
    end if

    mu6 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf))
    nu5 = ccradi(vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf))

    if ( max ( maxaft, nu5 ) <= max ( maxbef, mu6 ) * tolp1 ) then
      go to 50
    end if

    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  else

    if ( cc == a ) then
      call i4_swap ( a, b )
      call i4_swap ( va, vb )
    end if

    ind1 = htsrc(a,e,g,npt,sizht,fc,ht)

    if ( ind1 <= 0 ) then
      go to 80
    end if

    if ( fc(4,ind1) == b ) then
      h = fc(5,ind1)
    else
      h = fc(4,ind1)
    end if

    vh = vm(h)
    mu6 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,ve))
    mu7 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,vf))
    nu5 = ccradi(vcl(1,vg),vcl(1,vb),vcl(1,ve),vcl(1,vf))
    nu6 = ccradi(vcl(1,va),vcl(1,vh),vcl(1,ve),vcl(1,vf))
    nu7 = ccradi(vcl(1,vg),vcl(1,vh),vcl(1,ve),vcl(1,vf))

    if ( max ( maxaft, nu5, nu6, nu7 ) &
      <= max ( maxbef, mu6, mu7 ) * tolp1 ) then
      go to 50
    end if

    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0
    go to 60

  end if

50 continue

  if ( t1 .and. t2 ) then
    t1 = .false.
    first = .false.
    go to 30
  else
    go to 70
  end if

60 continue

  if ( msglvl == 4 ) then
    write ( *,630) 'combination swaps made'
  end if

  call swaptf ( top, npt, sizht, nfc, fc_max, vcl, vm, fc, ht, ntetra, &
    hdavfc, front, back, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  call swapes ( bndcon, 0, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
    ntetra, hdavfc, front, back, j, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  ind = ind + 1
  if ( nfc < ind ) then
    ind = 1
  end if

  ibeg = ind
  go to 10

70 continue

  ind = ind + 1

  if ( nfc < ind ) then
    ind = 1
  end if

  if ( ind /= ibeg ) then
    go to 10
  end if

  return

80 continue

  ierr = 300

  600 format (/1x,'imptrd')
  610 format (1x,'type ',a3,i7,' : ',5i7)
  620 format (4x,'face',3i7,' | ',2i7,' has type ',a3)
  630 format (4x,a)

  return
end
