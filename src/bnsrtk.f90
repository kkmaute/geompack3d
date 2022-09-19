subroutine bnsrtk ( k, binexp, n, a, map, bin, iwk, dx )

!*****************************************************************************80
!
!! BNSRTK bin sorts a set of KD points.
!
!  Discussion:
!
!    This routine uses bin sort to obtain the permutation of N K-dimensional
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
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, real ( kind = 8 ) BINEXP, the exponent for number of bins.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, the points of A with indices
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, the elements
!    are permuted so that
!      bin of MAP(1) <= bin of MAP(2) <= ... <= bin of MAP(N).
!
!    Workspace, integer BIN(1:N), used for bin numbers and permutation
!    of 1 to N.
!
!    Workspace, integer IWK(1:N), used for copy of MAP array.
!
!    Workspace, real ( kind = 8 ) DX(1:K), used for size of range of
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(k,*)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bin(n)
  real    ( kind = 8 ) binexp
  real    ( kind = 8 ) dx(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwk(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) nside
  integer ( kind = 4 ) nspow
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin

  nside = int(real(n)**(real(binexp)/real(k)) + 0.5D+00 )

  if ( nside <= 1 ) then
    return
  end if

  do i = 1, k
    xmin = minval ( a(i,map(1:n)) )
    xmax = maxval ( a(i,map(1:n)) )
    dx(i) = 1.0001D+00*(xmax - xmin)/ real ( nside, kind = 8 )
    if ( dx(i) == 0.0D+00 ) then
      dx(i) = 1.0D+00
    end if
  end do

  do i = 1, n

    j = map(i)
    iwk(i) = j
    map(i) = i
    b = int((a(1,j) - xmin) / dx(1))

    do l = 2, k

      m = int ( ( a(l,j) - xmin ) / dx(l) )

      if ( l == 2 ) then
        nspow = nside
      else
        nspow = nspow * nside
      end if

      if ( mod(m,2) == 0 ) then
        b = m*nspow + b
      else
        b = (m+1)*nspow - b - 1
      end if

    end do

    bin(i) = b

  end do

  call ihpsrt(1,n,1,bin,map)

  bin(1:n) = map(1:n)
  map(1:n) = iwk(bin(1:n))

  return
end
