subroutine mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

!*****************************************************************************80
!
!! MMASEP finds the best of four possible separators.
!
!  Discussion:
!
!    This routine finds the best of four possible separators according to
!    max-min angle criterion.
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
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter (in radians)
!    for accepting separator.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), coordinates of polygon
!    vertices in counterclockwise order where NVRT is number of vertices;
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, integer ( kind = 4 ) INDPVL(0:NVRT), indices in PVL of vertices;
!    INDPVL(I) = -K if ( XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, real ( kind = 8 ) IANG(1:*), interior angles.
!
!    Input, integer ( kind = 4 ) V(1:2), W(1:2), indices in XC, YC in range 0 to NVRT-1;
!    four possible separators are V(I),W(J), I,J = 1,2.
!
!    Output, integer ( kind = 4 ) I1, I2, indices in range 0 to NVRT-1 of best separator
!    according to max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) angle
  real    ( kind = 8 ) angmax
  real    ( kind = 8 ) angmin
  real    ( kind = 8 ) angtol
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) indpvl(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) v(2)
  integer ( kind = 4 ) w(2)
  real    ( kind = 8 ) xc(0:*)
  real    ( kind = 8 ) yc(0:*)

  tol = 100.0D+00 * epsilon ( tol )
  angmax = 0.0D+00

  do i = 1, 2

    l = v(i)
    k = indpvl(l)

    if ( 0 < k ) then
      alpha = iang(k)
    else
      alpha = pi
    end if

    do j = 1, 2

      m = w(j)

      if ( l == m ) then
        cycle
      end if

      k = indpvl(m)

      if ( 0 < k ) then
        beta = iang(k)
      else
        beta = pi
      end if

      gamma = angle ( xc(m), yc(m), xc(l), yc(l), xc(l+1), yc(l+1) )
      delta = angle ( xc(l), yc(l), xc(m), yc(m), xc(m+1), yc(m+1) )
      angmin = min ( gamma, alpha-gamma, delta, beta-delta )

      if ( angmax < angmin ) then
        angmax = angmin
        i1 = l
        i2 = m
      end if

    end do

  end do

  if ( angmax < angtol ) then
    i1 = -1
  end if

  return
end
