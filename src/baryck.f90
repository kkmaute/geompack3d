subroutine baryck ( k, ind, vcl, pt, alpha, degen, mat, ipvt )

!*****************************************************************************80
!
!! BARYCK computes the barycentric coordinates of a point in KD.
!
!  Discussion:
!
!    This routine computes the barycentric coordinates of a point with
!    respect to a simplex of K+1 vertices in K-dimensional space.
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
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, dimension of points and simplex,
!
!    Input, integer ( kind = 4 ) IND(1:K+1), indices in VCL of K-D vertices of simplex.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), K-D vertex coordinate list.
!
!    Input, real ( kind = 8 ) PT(1:K), K-D point for which barycentric
!    coordinates are to be computed.
!
!    Output, real ( kind = 8 ) ALPHA(1:K+1), barycentric coordinates
!    (if DEGEN = .FALSE.) such that
!    PT = ALPHA(1)*V[IND(1)] + ... + ALPHA(K+1)*V[IND(K+1)].
!
!    Output, logical DEGEN, is TRUE if the K+1 vertices form a
!    degenerate simplex.
!
!    Workspace, real ( kind = 8 ) MAT(1:K,1:K), matrix used for solving
!    system of linear equations.
!
!    Workspace, integer IPVT(1:K-1), pivot indices.
!
  implicit none

  integer ( kind = 4 ) k

  real    ( kind = 8 ) alpha(k+1)
  logical              degen
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k+1)
  integer ( kind = 4 ) ipvt(k-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mat(k,k)
  real    ( kind = 8 ) pt(k)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(k,*)

  tol = 100.0D+00 * epsilon ( tol )
  m = ind(k+1)

  do j = 1, k
    l = ind(j)
    do i = 1, k
      mat(i,j) = vcl(i,l) - vcl(i,m)
    end do
  end do

  alpha(1:k) = pt(1:k) - vcl(1:k,m)

  call lufac ( mat, k, k, tol, ipvt, degen )

  if ( .not. degen ) then
    call lusol ( mat, k, k, ipvt, alpha )
    alpha(k+1) = 1.0D+00 - sum ( alpha(1:k) )
  end if

  return
end
