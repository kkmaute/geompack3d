subroutine tmerge ( inter, nbl, ncr, chbl, chcr, ldv, vcl, til, tedg, &
  ierror )

!*****************************************************************************80
!
!! TMERGE forms triangles near the boundary by merging vertex chains.
!
!  Discussion:
!
!    This routine forms triangles in the strip near the boundary of a polygon
!    or inside the polygon by merging two chains of vertices.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, logical INTER, is TRUE iff at least one interior mesh vertex.
!
!    Input, integer ( kind = 4 ) NBL, the number of vertices on boundary cycle if INTER,
!    otherwise on left boundary chain.
!
!    Input, integer ( kind = 4 ) NCR, the number of vertices on closed walk if INTER,
!    otherwise on right boundary chain.
!
!    Input, integer ( kind = 4 ) CHBL(0:NBL), the indices in VCL of vertices on boundary
!    cycle or left boundary chain; if INTER, CHBL(NBL) = CHBL(0).
!
!    Input, integer ( kind = 4 ) CHCR(0:NCR), the indices in VCL of vertices on closed walk
!    or right boundary chain; if INTER, CHCR(NCR) = CHCR(0),
!    otherwise CHCR(0) is not referenced.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NT), the triangle incidence list, where NT =
!    NBL + NCR - K where K = 0 if INTER, else K = 2.
!
!    Output, integer ( kind = 4 ) TEDG(1:3,1:NT), the TEDG(J,I) refers to edge with vertices
!    TIL(J:J+1,I) and contains index of merge edge or NBL+NCR+1 for edge of
!    chains.  Note: It is assumed there is enough space in 2 arrays.
!
!    Output, integer ( kind = 4 ) IERROR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) nbl
  integer ( kind = 4 ) ncr

  integer ( kind = 4 ) chbl(0:nbl)
  integer ( kind = 4 ) chcr(0:ncr)
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibndry
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) in
  logical              inter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lri
  integer ( kind = 4 ) lrip1
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) tedg(3,nbl+ncr)
  integer ( kind = 4 ) til(3,nbl+ncr)
  real    ( kind = 8 ) vcl(ldv,*)
  real    ( kind = 8 ) xi
  real    ( kind = 8 ) xip1
  real    ( kind = 8 ) xj
  real    ( kind = 8 ) xjp1
  real    ( kind = 8 ) yi
  real    ( kind = 8 ) yip1
  real    ( kind = 8 ) yj
  real    ( kind = 8 ) yjp1

  ibndry = nbl + ncr + 1
  nt = 0

  if ( inter ) then

    nl = nbl
    nr = ncr
    i = 0
    j = 0

  else

    call mtredg ( .true., chbl(1), chcr(1), chbl(0), ibndry, nt, til, tedg )

    tedg(2,1) = ibndry
    if ( nbl + ncr <= 3 ) then
      return
    end if

    nl = nbl - 1
    nr = ncr - 1
    i = 1
    j = 1
    lri = 1
    lrip1 = 1

  end if
!
!  Main while loop for determining next triangle and edge.
!
  do

    if ( nl <= i .or. nr <= j ) then
      exit
    end if

    xi = vcl(1,chbl(i))
    yi = vcl(2,chbl(i))
    xip1 = vcl(1,chbl(i+1))
    yip1 = vcl(2,chbl(i+1))
    xj = vcl(1,chcr(j))
    yj = vcl(2,chcr(j))
    xjp1 = vcl(1,chcr(j+1))
    yjp1 = vcl(2,chcr(j+1))
    in = diaedg ( xjp1, yjp1, xj, yj, xi, yi, xip1, yip1 )

    if ( inter ) then
      lri = lrline ( xi, yi, xj, yj, xjp1, yjp1, 0.0D+00 )
      lrip1 = lrline ( xip1, yip1, xj, yj, xjp1, yjp1, 0.0D+00 )
    end if

    if ( in <= 0 .or. ( lri <= 0 .and. lrip1 <= 0 ) ) then

      call mtredg ( .true., chbl(i+1), chcr(j), chbl(i), ibndry, nt, til, tedg )
      i = i + 1

    else

      call mtredg ( .false., chbl(i), chcr(j+1), chcr(j), ibndry, nt, til, &
        tedg )
      j = j + 1

    end if

  end do
!
!  Add remaining triangles at end of strip or bottom of polygon.
!
  if ( i < nl ) then

    if ( ( .not. inter ) .and. j == nr ) then
      nl = nl + 1
    end if

    do

      call mtredg ( .true., chbl(i+1), chcr(j), chbl(i), ibndry, nt, til, tedg )
      i = i + 1

      if ( nl <= i ) then
        exit
      end if

    end do

  else
!
!  J < NR or I = NL = J = NR = 1
!
    if ( ( .not. inter ) .and. i == nl ) then
      nr = nr + 1
    end if

    do

      call mtredg ( .false., chbl(i), chcr(j+1), chcr(j), ibndry, nt, til, &
        tedg )

      if ( inter ) then

        lri = lrline ( vcl(1,chbl(i)), vcl(2,chbl(i)), &
          vcl(1,chcr(j+1)), vcl(2,chcr(j+1)), vcl(1,chcr(j)), &
          vcl(2,chcr(j)), 0.0D+00 )

        if ( 0 <= lri ) then
          ierror = 230
          return
        end if

      end if

      j = j + 1

      if ( nr <= j ) then
        exit
      end if

    end do

  end if

  if ( inter ) then
    if ( tedg(2,1) == 0 ) then
      tedg(2,1) = nbl + ncr
    else
      tedg(3,1) = nbl + ncr
    end if
  end if

  return
end
