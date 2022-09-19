subroutine bnsrt3 ( binexp, n, a, map, bin, iwk )

!*****************************************************************************80
!
!! BNSRT3 bin sorts a set of 3D points.
!
!  Discussion:
!
!    This routine uses bin sort to obtain the permutation of N 3-dimensional
!    double precision points so that points are in increasing bin
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
!    Input, real ( kind = 8 ) A(1:3,1:*), the array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, the points of A with indices
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, elements
!    are permuted so that:
!    bin of MAP(1) <= bin of MAP(2) <= ... <= bin of MAP(N)
!
!    Workspace, integer BIN(1:N), used for bin numbers and permutation
!    of 1 to N.
!
!    Workspace, integer IWK(1:N), used for copy of MAP array.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(3,*)
  integer ( kind = 4 ) bin(n)
  real    ( kind = 8 ) binexp
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  real    ( kind = 8 ) dz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwk(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) nside
  integer ( kind = 4 ) nsidsq
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) ymax
  real    ( kind = 8 ) ymin
  real    ( kind = 8 ) zmax
  real    ( kind = 8 ) zmin

  nside = int ( real(n)**( binexp / 3.0D+00 ) + 0.5D+00 )

  if ( nside <= 1 ) then
    return
  end if

  nsidsq = nside**2

  xmin = minval ( a(1,map(1:n)) )
  xmax = maxval ( a(1,map(1:n)) )
  ymin = minval ( a(2,map(1:n)) )
  ymax = minval ( a(2,map(1:n)) )
  zmin = minval ( a(3,map(1:n)) )
  zmax = minval ( a(3,map(1:n)) )

  dx = 1.0001D+00*(xmax - xmin)/ real ( nside, kind = 8 )
  dy = 1.0001D+00*(ymax - ymin)/ real ( nside, kind = 8 )
  dz = 1.0001D+00*(zmax - zmin)/ real ( nside, kind = 8 )
  if ( dx == 0.0D+00) dx = 1.0D+00
  if ( dy == 0.0D+00) dy = 1.0D+00
  if ( dz == 0.0D+00) dz = 1.0D+00

  do i = 1, n

    j = map(i)
    iwk(i) = j
    map(i) = i
    k = int((a(1,j) - xmin)/dx)
    l = int((a(2,j) - ymin)/dy)
    m = int((a(3,j) - zmin)/dz)

    if ( mod(l,2) == 0 ) then
      bin(i) = l*nside + m
    else
      bin(i) = (l+1)*nside - m - 1
    end if

    if ( mod(k,2) == 0 ) then
      bin(i) = k*nsidsq + bin(i)
    else
      bin(i) = (k+1)*nsidsq - bin(i) - 1
    end if

  end do

  call ihpsrt(1,n,1,bin,map)

  bin(1:n) = map(1:n)
  map(1:n) = iwk(bin(1:n))

  return
end
