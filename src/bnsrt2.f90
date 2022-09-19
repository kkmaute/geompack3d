subroutine bnsrt2 ( binexp, n, a, map, bin, iwk )

!*****************************************************************************80
!
!! BNSRT2 bin sorts a set of 2D points.
!
!  Discussion:
!
!    This routine uses a bin sort to obtain the permutation of N 2-dimensional
!    double precision points so that the points are in increasing bin
!    order, where the N points are assigned to about N**BINEXP bins.
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
!    Input, real ( kind = 8 ) BINEXP, the exponent for number of bins.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) A(1:2,1:*), array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, the points of A with
!    indices MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output,
!    elements are permuted so that bin of MAP(1) <= bin of MAP(2) <=
!    ... <= bin of MAP(N).
!
!    Workspace, integer BIN(1:N), used for bin numbers and permutation
!    of 1 to N.
!
!    Workspace, integer IWK(1:N), used for copy of MAP array.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(2,*)
  integer ( kind = 4 ) bin(n)
  real    ( kind = 8 ) binexp
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwk(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) nside
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) ymax
  real    ( kind = 8 ) ymin

  nside = int ( real ( n, kind = 8  )**( binexp / 2.0D+00 ) + 0.5D+00 )

  if ( nside <= 1 ) then
    return
  end if

  xmin = minval ( a(1,map(1:n)) )
  xmax = maxval ( a(1,map(1:n)) )
  ymin = minval ( a(2,map(1:n)) )
  ymax = minval ( a(2,map(1:n)) )

  dx = 1.0001D+00 * (xmax - xmin) / real ( nside, kind = 8 )
  dy = 1.0001D+00 * (ymax - ymin) / real ( nside, kind = 8 )

  if ( dx == 0.0D+00 ) then
    dx = 1.0D+00
  end if

  if ( dy == 0.0D+00 ) then
    dy = 1.0D+00
  end if

  do i = 1, n
    j = map(i)
    iwk(i) = j
    map(i) = i
    k = int((a(1,j) - xmin) / dx )
    l = int((a(2,j) - ymin) / dy )
    if ( mod(k,2) == 0 ) then
      bin(i) = k*nside + l
    else
      bin(i) = (k+1)*nside - l - 1
    end if
  end do

  call ihpsrt ( 1, n, 1, bin, map )

  bin(1:n) = map(1:n)
  map(1:n) = iwk(bin(1:n))

  return
end
