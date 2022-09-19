subroutine holvrt ( nhole, vcl, hvl, pvl, holv )

!*****************************************************************************80
!
!! HOLVRT determines top and bottom vertices of holes in polygonal regions.
!
!  Discussion:
!
!    This routine determines the top and bottom vertices of holes in
!    polygonal regions, and sorts the top vertices in decreasing (y,x) order
!    and bottom vertices in increasing (y,x) order.
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
!    Input, integer ( kind = 4 ) NHOLE, number of holes in the region.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NHOLE), head vertex list; HVL(I) is index in PVL of
!    head vertex of Ith hole.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), polygon vertex list; see routine DSPGDC.
!
!    Output, integer ( kind = 4 ) HOLV(1:NHOLE*2), indices in PVL of top and bottom
!    vertices of holes; first (last) NHOLE entries are for top (bottom)
!    vertices; top (bottom) vertices are sorted in decreasing
!    (increasing) lexicographic (y,x) order of coordinates.
!
  implicit none

  integer ( kind = 4 ) nhole

  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) holv(nhole*2)
  integer ( kind = 4 ) hv
  integer ( kind = 4 ) hvl(nhole)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) nhp1
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) y
  real    ( kind = 8 ) ymax
  real    ( kind = 8 ) ymin
!
!  Determine top and bottom vertices of holes.
!
  do i = 1, nhole

    hv = hvl(i)
    iv = hv

    do

      lv = pvl(loc,iv)

      if ( iv == hv ) then

        imin = iv
        imax = iv
        xmin = vcl(1,lv)
        ymin = vcl(2,lv)
        xmax = xmin
        ymax = ymin

      else

        x = vcl(1,lv)
        y = vcl(2,lv)

        if ( y < ymin .or. ( y == ymin .and. x < xmin ) ) then
          imin = iv
          xmin = x
          ymin = y
        else if ( ymax < y .or. ( y == ymax .and. xmax < x ) ) then
          imax = iv
          xmax = x
          ymax = y
        end if

      end if

      iv = pvl(succ,iv)

      if ( iv == hv ) then
        exit
      end if

    end do

    holv(i) = imax
    holv(i+nhole) = imin

  end do
!
!  Use linear insertion sort to sort the top vertices of holes
!  in decreasing (y,x) order, then bottom vertices in increasing
!  (y,x) order.  It is assumed NHOLE is small.
!
  do i = 2, nhole

    hv = holv(i)
    lv = pvl(loc,hv)
    x = vcl(1,lv)
    y = vcl(2,lv)
    j = i

30  continue

    iv = holv(j-1)
    lv = pvl(loc,iv)

    if ( vcl(2,lv) < y .or. ( y == vcl(2,lv) .and. vcl(1,lv) < x ) ) then
      holv(j) = iv
      j = j - 1
      if ( 1 < j ) then
        go to 30
      end if
    end if

    holv(j) = hv

  end do

  nhp1 = nhole + 1

  do i = nhp1+1, nhole+nhole

    hv = holv(i)
    lv = pvl(loc,hv)
    x = vcl(1,lv)
    y = vcl(2,lv)
    j = i

50  continue

    iv = holv(j-1)
    lv = pvl(loc,iv)

    if ( y < vcl(2,lv) .or. ( y == vcl(2,lv) .and. x < vcl(1,lv) ) ) then
      holv(j) = iv
      j = j - 1
      if ( nhp1 < j ) then
        go to 50
      end if
    end if

    holv(j) = hv

  end do

  return
end
