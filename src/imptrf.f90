subroutine imptrf ( bndcon, crit, npt, sizht, fc_max, vcl, vm, nfc, ntetra, &
  hdavfc, bf, fc,ht, ierr )

!*****************************************************************************80
!
!! IMPTRF improves a given triangulation in 3D.
!
!  Discussion:
!
!    This routine further improves a given 3D triangulation by applying
!    local transformations based on a local criterion.  Combination
!    swaps are used to remove poorly shaped tetrahedra.
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
!    Input, integer ( kind = 4 ) CRIT, the criterion code; 1 for (local max-min) solid angle
!    criterion, 2 for radius ratio criterion, 3 for mean ratio
!    criterion, 0 (or anything else) for no swaps.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, the size of the hash table HT; a good choice is a
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
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the  array of boundary face records;
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
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) crit
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) ff
  logical              first
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ierr
  logical              impr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) ind3
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) indy
  integer ( kind = 4 ) indz
  integer ( kind = 4 ) j
  real    ( kind = 8 ) minaft
  real    ( kind = 8 ) minbef
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mu1
  real    ( kind = 8 ) mu2
  real    ( kind = 8 ) mu3
  real    ( kind = 8 ) mu4
  real    ( kind = 8 ) mu5
  real    ( kind = 8 ) mu6
  real    ( kind = 8 ) mu7
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  real    ( kind = 8 ) nu1
  real    ( kind = 8 ) nu2
  real    ( kind = 8 ) nu3
  real    ( kind = 8 ) nu4
  real    ( kind = 8 ) nu5
  real    ( kind = 8 ) nu6
  real    ( kind = 8 ) nu7
  real    ( kind = 8 ) s(4)
  logical              t1
  logical              t2
  real    ( kind = 8 ) tetmu
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) top2
  character ( len = 3 ) typ2
  character ( len = 3 ) type
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
    write ( *,600) crit
  end if

  front = 0
  back = 0
  ibeg = 1
  ind = 1

10 continue

  if ( fc(1,ind) <= 0 .or. fc(5,ind) <= 0 ) then
    go to 70
  end if

  call ifacty ( ind, npt, sizht, vcl, vm, fc, ht, type, a, b, c, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( type /= 'N32' .and. type /= 'N44' ) then
    go to 70
  end if

  d = fc(4,ind)
  e = fc(5,ind)
  va = vm(a)
  vb = vm(b)
  vc = vm(c)
  vd = vm(d)
  ve = vm(e)

  if ( msglvl == 4 ) then
    write ( *,610) type,ind,a,b,c,d,e
  end if

  indx = htsrc ( a, b, d, npt, sizht, fc, ht )

  if ( indx <= 0 ) then
    go to 80
  end if

  if ( fc(4,indx) == c ) then
    f = fc(5,indx)
  else
    f = fc(4,indx)
  end if

  ind2 = htsrc(a,b,f,npt,sizht,fc,ht)

  if ( ind2 <= 0) go to 80

  if ( type == 'N44' ) then
    go to 20
  end if
!
!  TYPE == 'N32'
!
  if ( fc(4,ind2) /= e .and. fc(5,ind2) /= e ) then
    go to 70
  end if

  call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

  if ( msglvl == 4 ) then
    write ( *,620) aa,bb,cc,d,e,type
  end if

  if ( type /= 'T23' .and. type /= 'T32' .and. type /='T44' &
    .and. type /= 'N32' ) then
    go to 70
  end if

  vf = vm(f)
  mu1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
  mu2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
  mu3 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd),s)
  mu4 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,ve),s)
  nu1 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
  nu2 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)

  minbef = min ( mu1, mu2, mu3, mu4 )

  if ( type == 'T23' ) then

    minbef = minbef + tol

    if ( min ( nu1, nu2 ) <= minbef ) then
      go to 70
    end if

    nu3 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
    nu4 = tetmu(crit,vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)

    if ( max ( nu3, nu4 ) <= minbef ) then
      go to 70
    end if

    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0

    if ( min ( nu3, nu4 ) <= minbef ) then

      indy = htsrc(a,b,e,npt,sizht,fc,ht)

      if ( indy <= 0 ) then
        go to 80
      end if

      top2 = indx
      fc(7,indx) = indy
      fc(7,indy) = 0

      if ( nu4 < nu3 ) then
        call i4_swap ( a, b )
      end if

      call fndmsw(crit,npt,sizht,vcl,vm,fc,ht,a,b,d,e,f,minbef, &
        top,top2,impr,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      if ( .not. impr ) then
        go to 70
      end if

    end if

    go to 60

  else if ( type == 'T32' ) then

    if ( cc == a ) then
      call i4_swap ( va, vb )
    end if

    mu5 = tetmu(crit,vcl(1,vf),vcl(1,va),vcl(1,vd),vcl(1,ve),s)
    nu3 = tetmu(crit,vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve),s)

    if ( min ( nu1, nu2, nu3 ) <= min ( minbef, mu5 ) + tol ) then
      go to 70
    end if

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

    if ( ind1 <= 0 ) then
      go to 80
    end if

    if ( fc(4,ind1) == b ) then
      g = fc(5,ind1)
    else
      g = fc(4,ind1)
    end if

    vg = vm(g)
    mu5 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd),s)
    mu6 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve),s)
    nu3 = tetmu(crit,vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve),s)
    nu4 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
    nu5 = tetmu(crit,vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)

    if ( min ( nu1, nu2, nu3, nu4, nu5 ) <= &
      min ( minbef, mu5, mu6 ) + tol ) then
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

    ind3 = htsrc(a,f,d,npt,sizht,fc,ht)

    if ( ind3 <= 0 ) then
      go to 80
    end if

    if ( fc(4,ind3) == b ) then
      g = fc(5,ind3)
    else
      g = fc(4,ind3)
    end if

    ind1 = htsrc(a,f,g,npt,sizht,fc,ht)

    if ( ind1 <= 0 ) then
      go to 80
    end if

    if ( fc(4,ind1) /= e .and. fc(5,ind1) /= e ) then
      go to 70
    end if

    call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

    if ( type /= 'T23' ) then
      go to 70
    end if

    vg = vm(g)
    mu5 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd),s)
    mu6 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve),s)
    minbef = min(minbef,mu5,mu6) + tol
    nu3 = tetmu(crit,vcl(1,vb),vcl(1,vd),vcl(1,ve),vcl(1,vf),s)
    if ( min(nu1,nu2,nu3) <= minbef) go to 70
    nu4 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
    nu5 = tetmu(crit,vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
    if ( max ( nu4, nu5 ) <= minbef ) go to 70
    top = ind1
    fc(7,ind1) = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0

    if ( min ( nu4, nu5 ) <= minbef ) then

      indy = htsrc(a,b,e,npt,sizht,fc,ht)
      indz = htsrc(a,f,e,npt,sizht,fc,ht)
      if ( indy <= 0 .or. indz <= 0) go to 80
      top2 = indx
      fc(7,indx) = indy
      fc(7,indy) = indz
      fc(7,indz) = ind3
      fc(7,ind3) = 0

      if ( nu5 < nu4 ) then
        call i4_swap ( a, f )
      end if

      call fndmsw(crit,npt,sizht,vcl,vm,fc,ht,a,f,d,e,g,minbef, &
        top,top2,impr,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      if ( .not. impr ) then
        go to 70
      end if

    end if

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

  if ( ind2 <= 0 ) then
    go to 80
  end if

  if ( fc(4,ind2) /= e .and. fc(5,ind2) /= e ) then
    go to 70
  end if

  call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

  if ( msglvl == 4 ) then
    write ( *,620) aa,bb,cc,e,f,type
  end if

  t1 = (type == 'T23' .or. type == 'T32' .or. type =='T44')
  call ifacty(ind3,npt,sizht,vcl,vm,fc,ht,typ2,aa,bb,c2,ierr)

  if ( msglvl == 4 ) then
    write ( *,620) aa,bb,c2,d,g,typ2
  end if

  t2 = (typ2 == 'T23' .or. typ2 == 'T32' .or. typ2 =='T44')
  if ( ( .not. t1 ) .and. ( .not. t2 ) ) go to 70
  first = .true.

30 continue

  if ( .not. t1 ) then

    j = ind2
    ind2 = ind3
    ind3 = j
    type = typ2
    cc = c2

    call i4_swap ( d, e )
    call i4_swap ( vd, ve )
    call i4_swap ( f, g )

  end if

40 continue

  vf = vm(f)
  vg = vm(g)

  if ( first ) then
    mu1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
    mu2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
    mu3 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd),s)
    mu4 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,ve),s)
    mu5 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,vf),s)
    nu1 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
    nu2 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
    nu3 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
    nu4 = tetmu(crit,vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
  else
    nu3 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
    nu4 = tetmu(crit,vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
  end if

  if ( first ) then
    minbef = min(mu1,mu2,mu3,mu4,mu5)
    minaft = min(nu1,nu2,nu3,nu4)
  else
    minaft = min(nu1,nu2,nu3,nu4)
  end if

  if ( type == 'T23' ) then

    if ( minaft <= minbef + tol) go to 50
    nu5 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
    nu6 = tetmu(crit,vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
    if ( max ( nu5, nu6 ) <= minbef + tol ) go to 50
    top = ind2
    fc(7,ind2) = ind
    fc(7,ind) = 0

    if ( min ( nu5, nu6 ) <= minbef + tol ) then

      indy = htsrc(a,b,e,npt,sizht,fc,ht)
      if ( indy <= 0) go to 80
      top2 = ind3
      fc(7,ind3) = indy
      fc(7,indy) = 0

      if ( nu5 <= nu6 ) then
        aa = a
        bb = b
      else
        aa = b
        bb = a
      end if

      dd = e
      ee = f
      ff = g
      call fndmsw(crit,npt,sizht,vcl,vm,fc,ht,aa,bb,dd,ee,ff, &
        minbef+tol,top,top2,impr,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      if ( .not. impr) go to 50

    end if

    go to 60

  else if ( type == 'T32' ) then

    if ( cc == a ) then
      call i4_swap ( a, b )
      call i4_swap ( va, vb )
    end if

    mu6 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
    nu5 = tetmu(crit,vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
    if ( min(minaft,nu5) <= min(minbef,mu6) + tol) go to 50
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
    if ( ind1 <= 0) go to 80

    if ( fc(4,ind1) == b ) then
      h = fc(5,ind1)
    else
      h = fc(4,ind1)
    end if

    vh = vm(h)
    mu6 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,ve),s)
    mu7 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,vf),s)
    nu5 = tetmu(crit,vcl(1,vg),vcl(1,vb),vcl(1,ve),vcl(1,vf),s)
    nu6 = tetmu(crit,vcl(1,va),vcl(1,vh),vcl(1,ve),vcl(1,vf),s)
    nu7 = tetmu(crit,vcl(1,vg),vcl(1,vh),vcl(1,ve),vcl(1,vf),s)

    if ( min ( minaft, nu5, nu6, nu7 ) <= min ( minbef, mu6, mu7 ) + tol ) then
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

  call swaptf(top,npt,sizht,nfc,fc_max,vcl,vm,fc,ht,ntetra,hdavfc, &
    front,back, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  call swapmu(bndcon,crit,npt,sizht,nfc,fc_max,vcl,vm,bf,fc,ht, &
    ntetra,hdavfc,front,back,j, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  ind = ind + 1

  if ( nfc < ind ) then
    ind = 1
  end if

  ibeg = ind
  go to 10

70  continue

  ind = ind + 1
  if ( nfc < ind ) then
    ind = 1
  end if

  if ( ind /= ibeg) go to 10
  return

80 continue

  ierr = 300

  600 format (/1x,'imptrf: criterion =',i3)
  610 format (1x,'type ',a3,i7,' : ',5i7)
  620 format (4x,'face',3i7,' | ',2i7,' has type ',a3)
  630 format (4x,a)

  return
end
