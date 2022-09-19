subroutine fndsph ( xh, yh, nvrt, xc, yc, ivis, theta, nv, iv, x, y, link, &
  angsep, v )

!*****************************************************************************80
!
!! FNDSPH finds a separator from top or bottom hole vertex.
!
!  Discussion:
!
!    This routine finds a separator from top or bottom hole vertex (XH,YH)
!    using a max-min angle criterion from list of vertices in
!    increasing polar angle with respect to horizontal ray through (XH,YH).
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
!    Input, real ( kind = 8 ) XH, YH, the coordinates of hole vertex.
!
!    Input, integer ( kind = 4 ) NVRT, (number of vertices) - 1.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates
!    of possible endpoints of a separator.
!
!    Input, integer ( kind = 4 ) IVIS(0:NVRT), contains information about the vertices of
!    XC, YC arrays with respect to X, Y arrays; IVIS(I) = J greater than 0 if
!    (XC(I),YC(I)) = (X(K),Y(K)) where K = LINK(J), and
!    IVIS(I) = J < 0 if ( XC(I),YC(I)) lies in interior of
!    edge starting at (X(-J),Y(-J)).
!
!    Input, real ( kind = 8 ) THETA(0:NVRT), the polar angles of vertices
!    in increasing order with respect to horizontal ray through (XH,YH);
!    THETA(NVRT) = PI.
!
!    Input, integer ( kind = 4 ) NV, (number of vertices to be considered as endpoint of a
!    separator) - 1.
!
!    Input, integer ( kind = 4 ) IV(0:NV), the indices in increasing order of vertices
!    in XC, YC arrays to be considered as separator endpoint; angle
!    between consecutive vertices is assumed to be < 180
!    degrees; it is also assumed that 0, NVRT not in array.
!
!    Input, X(1:*), Y(1:*), LINK(1:*), used for 2D representation of
!    decomposition of multiply-connected polygonal face
!    of hole polygons; see routine SPDECH.
!
!    Output, real ( kind = 8 ) ANGSEP, the minimum of 4 angles at boundary
!    resulting from separator.
!
!    Output, integer ( kind = 4 ) V, the index of separator endpoint in XC, YC arrays.
!
  implicit none

  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) angle
  real    ( kind = 8 ) angmin
  real    ( kind = 8 ) angsep
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iv(0:nv)
  integer ( kind = 4 ) ivis(0:nvrt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) link(*)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta(0:nvrt)
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) v
  real    ( kind = 8 ) x(*)
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xh
  real    ( kind = 8 ) y(*)
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) yh

  angsep = 0.0D+00
  tol = 100.0D+00 * epsilon ( tol )

  do k = 0, nv

    i = iv(k)
    j = ivis(i)

    if ( 0 < j ) then
      l = link(link(j))
      alpha = angle(x(j),y(j),xc(i),yc(i),xh,yh)
      beta = angle(xh,yh,xc(i),yc(i),x(l),y(l))
    else
      j = -j
      alpha = angle(x(j),y(j),xc(i),yc(i),xh,yh)
      beta = pi - alpha
    end if

    angmin = min ( theta(i), pi - theta(i), alpha, beta )

    if ( angsep < angmin ) then
      angsep = angmin
      v = i
    end if

  end do

  return
end
