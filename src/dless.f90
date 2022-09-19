function dless ( k, p, q )

!*****************************************************************************80
!
!! DLESS determines the lexicographically lesser of two double precision values.
!
!  Discussion:
!
!    This routine determines whether P is lexicographically less than Q in
!    floating point arithmetic.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, dimension of points.
!
!    Input, real ( kind = 8 ) P(1:K), Q(1:K), two points.
!
!    Output, logical DLESS, TRUE if P < Q, FALSE otherwise.
!
  implicit none

  real    ( kind = 8 ) cmax
  logical              dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real    ( kind = 8 ) p(k)
  real    ( kind = 8 ) q(k)
  real    ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( abs ( p(i) - q(i) ) <= tol * cmax .or. cmax <= tol ) then
      cycle
    end if

     if ( p(i) < q(i) ) then
       dless = .true.
     else
       dless = .false.
     end if

     return

  end do

  dless = .false.

  return
end
