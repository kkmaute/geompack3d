subroutine ccsphk ( k, intest, ind, vcl, pt, center, radsq, in, mat, ipvt )

!*****************************************************************************80
!
!! CCSPHK finds the circumsphere through a simplex in KD.
!
!  Discussion:
!
!    This routine finds the center and the square of the radius of the
!    circumsphere through K+1 vertices of a K-D simplex, and possibly
!    determines whether another K-D point is inside the (hyper)sphere.
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
!    Input, integer ( kind = 4 ) K, the dimension of points and simplex.
!
!    Input, logical INTEST, TRUE iff test for point PT in sphere to be made.
!
!    Input, integer ( kind = 4 ) IND(1:K+1), the indices in VCL of K-D vertices of simplex.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), K-D vertex coordinate list.
!
!    Input, real ( kind = 8 ) PT(1:K), the point for which sphere test
!    is applied (referenced iff INTEST is TRUE).
!
!    Output, real ( kind = 8 ) CENTER(1:K), the center of circumsphere;
!    undefined if K+1 vertices of simplex lie on same K-D hyperplane.
!
!    Output, real ( kind = 8 ) RADSQ, the square of radius of sphere;
!    -1 in degenerate case.
!
!    Output, integer ( kind = 4 ) IN, contains following value if INTEST is TRUE:
!     2 if degenerate simplex;
!     1 if PT inside sphere;
!     0 if PT on sphere;
!    -1 if PT outside sphere
!
!    Workspace, real ( kind = 8 ) MAT(1:K,1:K), matrix used for solving
!    system of linear equations.
!
!    Workspace, integer IPVT(1:K-1), pivot indices
!
  implicit none

  integer ( kind = 4 ) k

  real    ( kind = 8 ) center(k)
  real    ( kind = 8 ) dsq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind(k+1)
  logical              intest
  integer ( kind = 4 ) ipvt(k-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mat(k,k)
  real    ( kind = 8 ) pt(k)
  real    ( kind = 8 ) radsq
  logical              singlr
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(k,*)

  tol = 100.0D+00 * epsilon ( tol )
  m = ind(k+1)

  do i = 1, k
    l = ind(i)
    sum2 = 0.0D+00
    do j = 1, k
      mat(i,j) = vcl(j,l) - vcl(j,m)
      sum2 = sum2 + mat(i,j)**2
    end do
    center(i) = 0.5D+00*sum2
  end do

  call lufac ( mat, k, k, tol, ipvt, singlr )

  if ( singlr ) then
    in = 2
    radsq = -1.0D+00
  else
    call lusol(mat,k,k,ipvt,center)
    radsq = 0.0D+00
    do i = 1, k
      radsq = radsq + center(i)**2
      center(i) = center(i) + vcl(i,m)
    end do
  end if

  if ( intest .and. .not. singlr ) then

   dsq = 0.0D+00
   do i = 1, k
     dsq = dsq + (pt(i) - center(i))**2
   end do

   if ( ( 1.0D+00 + tol ) * radsq < dsq ) then
     in = -1
   else if ( dsq < (1.0D+00 - tol)*radsq ) then
     in = 1
   else
     in = 0
   end if

  end if

  return
end
