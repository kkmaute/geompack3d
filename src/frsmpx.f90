subroutine frsmpx ( k, shift, nv, vcl, map, inds, ipvt, mat, ierr )

!*****************************************************************************80
!
!! FRSMPX shifts vertices to the first K+1 are in general position in KD.
!
!  Discussion:
!
!    This routine shifts or swaps vertices if necessary so that the
!    first K+1 vertices are not in same hyperplane (so first simplex is valid).
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
!    Input, integer ( kind = 4 ) K, dimension of points.
!
!    Input, logical SHIFT, if TRUE, MAP(3:K+1) may be updated due to shift,
!    else they may be updated due to swaps; in former case,
!    it is assumed MAP gives vertices in lexicographic order.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) MAP(1:NV).  On input, contains vertex indices of VCL.
!    On output, shifted or K-1 swaps applied if necessary so vertices
!    indexed by MAP(1), ..., MAP(K+1) not in same hyperplane.
!
!    Output, integer ( kind = 4 ) INDS(3:K+1), indices such that
!    MAP_in(INDS(I)) = MAP_out(I).
!
!    Workspace, real ( kind = 8 ) MAT(1:K,1:K), matrix used for
!    determining rank.
!
!    Workspace, integer IPVT(1:K), pivot indices.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) c
  real    ( kind = 8 ) cmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inds(3:k+1)
  integer ( kind = 4 ) ipvt(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) map(nv)
  real    ( kind = 8 ) mat(k,k)
  real    ( kind = 8 ) mult
  real    ( kind = 8 ) pivot
  real    ( kind = 8 ) rtol
  logical              shift
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(k,*)
!
!  First check that consecutive vertices are not identical.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( shift ) then
    l = nv - 1
  else
    l = 1
  end if

  m1 = map(1)

  do i = 1, l

    m = m1
    m1 = map(i+1)

    do j = 1, k
      cmax = max ( abs ( vcl(j,m) ), abs ( vcl(j,m1) ) )
      if ( tol * cmax < abs ( vcl(j,m) - vcl(j,m1) ) .and. tol < cmax ) then
        go to 20
      end if
    end do

    ierr = 402
    return

20  continue

  end do
!
!  Find indices INDS(3), ..., INDS(K+1).
!
  m1 = map(1)
  cmax = 0.0D+00
  c = 1
  i = 1

30 continue

  i = i + 1

  if ( nv < i ) then
    ierr = 403
    return
  end if

  m = map(i)

  do j = 1, k
    mat(j,c) = vcl(j,m) - vcl(j,m1)
    cmax = max ( cmax, abs ( mat(j,c) ) )
  end do

  rtol = tol*cmax

  do jj = 1, c-1
    l = ipvt(jj)
    mult = mat(l,c)
    mat(l,c) = mat(jj,c)
    do ii = jj+1, k
      mat(ii,c) = mat(ii,c) - mult * mat(ii,jj)
    end do
  end do

  l = c

  do j = c+1, k
    if ( abs ( mat(l,c) ) < abs ( mat(j,c) ) ) then
      l = j
    end if
  end do

  pivot = mat(l,c)

  if ( 1 < c ) then
    if ( abs ( pivot ) < rtol ) go to 30
    inds(c+1) = i
  end if

  ipvt(c) = l

  if ( l /= c ) then
    mat(l,c) = mat(c,c)
    mat(c,c) = pivot
  end if

  do ii = c+1, k
    mat(ii,c) = mat(ii,c) / pivot
  end do

  if ( c < k ) then
    c = c + 1
    go to 30
  end if
!
!  Shift or swap elements of MAP if necessary.
!
  if ( shift ) then

    do i = 3, k+1
      if ( i < inds(i) ) then
        ipvt(i-2) = map(inds(i))
      end if
    end do

    do i = k+1, 3, -1

      if ( i < inds(i) ) then

        l = k + 2 - i

        if ( 3 < i ) then
          m = inds(i-1) + 1
        else
          m = 3
        end if

        do j = inds(i)-1, m, -1
          map(j+l) = map(j)
        end do

      end if

    end do

    do i = 3, k+1
      if ( i < inds(i) ) then
        map(i) = ipvt(i-2)
      end if
    end do

  else

    do i = 3, k+1

      if ( i < inds(i) ) then
        m = map(i)
        map(i) = map(inds(i))
        map(inds(i)) = m
      end if

    end do

  end if

  return
end
