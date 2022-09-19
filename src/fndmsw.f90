subroutine fndmsw ( crit, npt, sizht, vcl, vm, fc, ht, a, b, d, e, f, minbef, &
  top, top2, impr, ierr )

!*****************************************************************************80
!
!! FNDMSW finds local transformation that improve a 3D triangulation.
!
!  Discussion:
!
!    This routine finds a sequence of 3 or more local transformations to
!    improve a 3D triangulation.
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
!    Input, integer ( kind = 4 ) CRIT, criterion code; 1 for (local max-min) solid angle
!    criterion, 2 for radius ratio criterion, 3 for mean ratio
!    criterion, 0 (or anything else) for no swaps.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine
!    DTRIS3.  On output, some faces may be added to the list.
!
!    Input, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) A,B,D,E,F, the vertices (local indices) in
!    configuration with T23 face AFB|DE swapped out to produce only 1
!    tetra AFDE with measure <= MINBEF; try to apply a T32 swap to remove AFDE
!    where AF is the desired edge to be removed.
!
!    Input, real ( kind = 8 ) MINBEF, the (min tetra measure of an existing
!    tetra of swap) + TOL.
!
!    Input/output, integer ( kind = 4 ) TOP; on input, a pointer to a list of 2 or 3
!    faces to be possibly swapped.  On output, if IMPR is TRUE, TOP is
!    a pointer to a list of faces to be swapped, but if IMPR is FALSE,
!    TOP is 0 and all faces were removed from the list.
!
!    Input/output, integer ( kind = 4 ) TOP2.  On input, pointer to stack of other
!    faces of T32 or T44 swaps.  On output, is set to zero, and stack
!    is emptied.
!
!    Output, logical IMPR, is TRUE iff improvement is possible.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  real    ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) crit
  integer ( kind = 4 ) d
  logical              degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  logical              impr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) indy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  integer ( kind = 4 ), parameter :: maxtf = 13
  real    ( kind = 8 ) minbef
  real    ( kind = 8 ) nu1
  real    ( kind = 8 ) nu2
  real    ( kind = 8 ) nu3
  real    ( kind = 8 ) nu4
  real    ( kind = 8 ) nu5
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) tetmu
  integer ( kind = 4 ) top
  integer ( kind = 4 ) top2
  character ( len = 3 ) typ2
  character ( len = 3 ) type
  integer ( kind = 4 ) va
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vf
  integer ( kind = 4 ) vg
  integer ( kind = 4 ) vh
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  kf = 2
  impr = .false.
  ptr = top

10 continue

  if ( maxtf <= kf ) then
    go to 40
  end if

  indx = htsrc ( a, f, d, npt, sizht, fc, ht )

  if ( indx <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(4,indx) == b ) then
    g = fc(5,indx)
  else
    g = fc(4,indx)
  end if

  ind = htsrc ( a, f, e, npt, sizht, fc, ht )

  if ( ind <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(4,ind) == b ) then
    h = fc(5,ind)
  else
    h = fc(4,ind)
  end if

  if ( g <= 0 .or. h <= 0 ) then
    go to 40
  end if

  if ( g /= h ) then

    ind1 = htsrc ( a, f, h, npt, sizht, fc, ht )

    if ( ind1 <= 0 ) then
      ierr = 300
      return
    end if

    if ( fc(4,ind1) /= g .and. fc(5,ind1) /= g ) then
      go to 40
    end if

    call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

    if ( type /= 'T23' ) then

      ind2 = ind1
      c2 = cc
      typ2 = type
      ind1 = htsrc ( a, f, g, npt, sizht, fc, ht )

      if ( ind1 <= 0 ) then
        ierr = 300
        return
      end if

      call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

      if ( type == 'T23' .or. type == 'N32' ) then
        call i4_swap ( d, e )
        call i4_swap ( g, h )
      else if ( typ2 == 'N32' ) then
        type = typ2
        ind1 = ind2
        cc = c2
      else
        go to 40
      end if

    end if

  end if

  va = vm(a)
  vd = vm(d)
  ve = vm(e)
  vf = vm(f)
  vg = vm(g)

  call baryth(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),vcl(1,vg), &
    alpha,degen)

  if ( degen ) then
    ierr = 301
    return
  end if

  if ( 0.0D+00 < alpha(4) ) then
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

  if ( kneg /= 2 .or. kzero /= 0 .or. 0.0D+00 <= alpha(3) ) then
    go to 40
  end if

  if ( fc(7,ind) /= -1 .or. fc(7,indx) /= -1) then
    go to 40
  end if

  indy = htsrc(a,f,g,npt,sizht,fc,ht)

  if ( indy <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(7,indy) /= -1 ) then
    go to 40
  end if

  nu1 = tetmu(crit,vcl(1,va),vcl(1,vd),vcl(1,ve),vcl(1,vg),alpha)
  nu2 = tetmu(crit,vcl(1,vf),vcl(1,vd),vcl(1,ve),vcl(1,vg),alpha)

  if ( min ( nu1, nu2 ) <= minbef ) then
    go to 40
  end if

  fc(7,ind) = fc(7,ptr)
  fc(7,ptr) = ind
  kf = kf + 1
!
!  Last face added to middle of list has type T32.
!
  if ( g == h ) then
    impr = .true.
    go to 50
  end if

  fc(7,indx) = indy
  fc(7,indy) = top2
  top2 = indx

  if ( type == 'N32') then
    go to 30
  end if

  if ( fc(7,ind1) /= -1) then
    go to 40
  end if

  vh = vm(h)
  nu3 = tetmu(crit,vcl(1,va),vcl(1,vh),vcl(1,ve),vcl(1,vg),alpha)
  nu4 = tetmu(crit,vcl(1,vf),vcl(1,vh),vcl(1,ve),vcl(1,vg),alpha)

  if ( max ( nu3, nu4 ) <= minbef ) then
    go to 40
  end if

  fc(7,ind1) = fc(7,ptr)
  fc(7,ptr) = ind1
  ptr = ind1
  kf = kf + 1
!
!  Last face added to middle of list has type T23.
!
  if ( minbef < min ( nu3, nu4 ) ) then
    impr = .true.
    go to 50
  end if

  if ( nu4 < nu3 ) then
    call i4_swap ( a, f )
  end if

  b = f
  d = g
  f = h
  go to 10

30 continue

  if ( cc == a ) then
    call i4_swap ( a, f )
    call i4_swap ( va, vf )
  end if

  indx = htsrc(a,h,e,npt,sizht,fc,ht)

  if ( indx <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(7,indx) /= -1 ) then
    go to 40
  end if

  if ( fc(4,indx) == f ) then
    i = fc(5,indx)
  else
    i = fc(4,indx)
  end if

  ind2 = htsrc(a,h,i,npt,sizht,fc,ht)

  if ( ind2 <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(4,ind2) /= g .and. fc(5,ind2) /= g ) then
    go to 40
  end if

  call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc,ierr)

  if ( type /= 'T23' ) then
    go to 40
  end if

  vh = vm(h)
  vi = vm(i)
  nu3 = tetmu(crit,vcl(1,ve),vcl(1,vf),vcl(1,vg),vcl(1,vh),alpha)

  if ( nu3 <= minbef) go to 40

  nu4 = tetmu(crit,vcl(1,va),vcl(1,vi),vcl(1,ve),vcl(1,vg),alpha)
  nu5 = tetmu(crit,vcl(1,vh),vcl(1,vi),vcl(1,ve),vcl(1,vg),alpha)

  if ( max ( nu4, nu5 ) <= minbef ) then
    go to 40
  end if

  if ( fc(7,ind1) /= -1 .or. fc(7,ind2) /= -1) go to 40

  indy = htsrc(a,h,g,npt,sizht,fc,ht)

  if ( indy <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(7,indy) /= -1) go to 40

  fc(7,ind1) = fc(7,ptr)
  fc(7,ptr) = ind2
  fc(7,ind2) = ind1
  ptr = ind2
  kf = kf + 2
!
!  Last 2 faces added to middle of list have type T23, T32.
!
  if ( minbef < min ( nu4, nu5 ) ) then
    impr = .true.
    go to 50
  end if

  fc(7,indx) = indy
  fc(7,indy) = top2
  top2 = indx

  if ( nu5 < nu4 ) then
    call i4_swap ( a, h )
  end if

  b = h
  d = g
  f = i
  go to 10

40 continue

  do

    ptr = top
    top = fc(7,ptr)
    fc(7,ptr) = -1

    if ( top == 0 ) then
      exit
    end if

  end do

50 continue

  do

    ptr = top2
    top2 = fc(7,ptr)
    fc(7,ptr) = -1

    if ( top2 == 0 ) then
      exit
    end if

  end do

  return
end
