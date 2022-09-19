function opsidk ( k, ind, vcl, eflag, pta, ptb, mat, vec )

!*****************************************************************************80
!
!! OPSIDK tests if points are on opposite sides of a face in KD.
!
!  Discussion:
!
!    This routine tests if points PTA, PTB are on opposite sides of (K-1)-D
!    face formed by K K-D vertices.
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
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) IND(1:K), the indices in VCL of K-D vertices of face.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), the K-D vertex coordinate list.
!
!    Input, logical EFLAG, TRUE if PTA = PTB: used to determine whether
!    PTA lies in hyperplane containing face, only PTA referenced.
!
!    Input, real ( kind = 8 ) PTA(1:K), PTB(1:K), the K-D points for which
!    test is applied.
!
!    Workspace, real ( kind = 8 ) MAT(1:K-1,1:K), the matrix used for
!    solving system of homogeneous linear equations.
!
!    Workspace, real ( kind = 8 ) VEC(1:K), the vector used for
!    hyperplane normal.
!
!    Output, integer ( kind = 4 ) OPSIDK,
!    +1 if PTA, PTB on opposite sides;
!    -1 if on same side;
!     0 if face is degenerate, or PTA or PTB is on same hyperplane as face.
!
  implicit none

  integer ( kind = 4 ) k

  real    ( kind = 8 ) dpa
  real    ( kind = 8 ) dpb
  logical              eflag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mat(k-1,k)
  integer ( kind = 4 ) opsidk
  real    ( kind = 8 ) pta(k)
  real    ( kind = 8 ) ptb(k)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) vcl(k,*)
  real    ( kind = 8 ) vec(k)

  tol = 100.0D+00 * epsilon ( tol )
  opsidk = 0
  km1 = k - 1
  m = ind(k)

  t = 0.0D+00
  do i = 1, km1
    l = ind(i)
    do j = 1, k
      mat(i,j) = vcl(j,l) - vcl(j,m)
      t = max ( t, abs ( mat(i,j) ) )
    end do
  end do

  tolabs = tol * t
!
!  Use Gaussian elimination with partial pivoting to solve K-1 by K
!  homogeneous system. Face is degenerate if 2 zero pivots occur.
!
  r = k
  s = 0

  do l = 1, k-2

30  continue

    ll = l + s
    m = l
    do i = l+1, km1
      if ( abs ( mat(m,ll) ) < abs ( mat(i,ll) ) ) then
        m = i
      end if
    end do

    t = mat(m,ll)
    mat(m,ll) = mat(l,ll)
    mat(l,ll) = t

    if ( abs ( t ) <= tolabs ) then
      if ( s == 1 ) then
        return
      else
        r = l
        s = 1
        go to 30
      end if
    end if

    do i = l+1, km1
      mat(i,ll) = mat(i,ll) / t
    end do

    do j = ll+1, k
      t = mat(m,j)
      mat(m,j) = mat(l,j)
      mat(l,j) = t
      do i = l+1, km1
        mat(i,j) = mat(i,j) - mat(i,ll)*t
      end do
    end do

  end do

  if ( abs ( mat(km1,km1+s) ) <= tolabs ) then

    if ( s == 1 ) then
      return
    else
      if ( abs ( mat(km1,k) ) <= tolabs ) then
        return
      end if
      r = km1
    end if

  end if
!
!  Matrix has full rank. If R <= K-1 then column R has a zero pivot
!  and VEC(R+1:K) = 0. Use VEC(1:R) for opposite test.
!
  vec(r) = -1.0D+00
  do l = r-1, 1, -1
    t = mat(l,r) / mat(l,l)
    vec(l) = t
    do i = 1, l-1
      mat(i,r) = mat(i,r) - mat(i,l)*t
    end do
  end do

  m = ind(k)
  dpa = 0.0D+00

  if ( eflag ) then

    do i = 1, r
      dpa = dpa + vec(i)*(pta(i) - vcl(i,m))
    end do

    if ( abs ( dpa ) <= tolabs ) then
      opsidk = 0
    else
      opsidk = -1
    end if

  else

    dpb = 0.0D+00
    do i = 1, r
      dpa = dpa + vec(i) * ( pta(i) - vcl(i,m) )
      dpb = dpb + vec(i) * ( ptb(i) - vcl(i,m) )
    end do

    if ( abs ( dpa ) <= tolabs .or. abs ( dpb ) <= tolabs ) then
      opsidk = 0
    else if ( dpa * dpb < 0.0D+00 ) then
      opsidk = 1
    else
      opsidk = -1
    end if

  end if

  return
end
