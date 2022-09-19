subroutine insph ( a, b, c, d, center, rad )

!*****************************************************************************80
!
!! INSPH finds the center and radius of the insphere of a tetrahedron.
!
!  Discussion:
!
!    This routine finds the center and radius of the insphere of a tetrahedron.
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
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), 4 vertices of
!    tetrahedron.
!
!    Output, real ( kind = 8 ) CENTER(1:3), the center of insphere; undefined
!    if A,B,C,D coplanar.
!
!    Output, real ( kind = 8 ) RAD, the radius of insphere; 0 if A,B,C,D
!    coplanar.
!
  implicit none

  real    ( kind = 8 ) a(3)
  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) ad(3)
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) bc(3)
  real    ( kind = 8 ) bd(3)
  real    ( kind = 8 ) c(3)
  real    ( kind = 8 ) center(3)
  real    ( kind = 8 ) d(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) mat(4,4)
  integer ( kind = 4 ) r
  real    ( kind = 8 ) rad
  real    ( kind = 8 ) rtol
  logical              singlr
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
!
!  Compute unit outward (or inward) normals and equations of 4 faces.
!
  tol = 100.0D+00 * epsilon ( tol )
  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)
  ad(1:3) = d(1:3) - a(1:3)
  bc(1:3) = c(1:3) - b(1:3)
  bd(1:3) = d(1:3) - b(1:3)

  rtol = tol * max ( abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
                     abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
                     abs ( c(1) ), abs ( c(2) ), abs ( c(3) ), &
                     abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  mat(1,1) = ac(2) * ab(3) - ac(3) * ab(2)
  mat(1,2) = ac(3) * ab(1) - ac(1) * ab(3)
  mat(1,3) = ac(1) * ab(2) - ac(2) * ab(1)
  mat(2,1) = ab(2) * ad(3) - ab(3) * ad(2)
  mat(2,2) = ab(3) * ad(1) - ab(1) * ad(3)
  mat(2,3) = ab(1) * ad(2) - ab(2) * ad(1)
  mat(3,1) = ad(2) * ac(3) - ad(3) * ac(2)
  mat(3,2) = ad(3) * ac(1) - ad(1) * ac(3)
  mat(3,3) = ad(1) * ac(2) - ad(2) * ac(1)
  mat(4,1) = bc(2) * bd(3) - bc(3) * bd(2)
  mat(4,2) = bc(3) * bd(1) - bc(1) * bd(3)
  mat(4,3) = bc(1) * bd(2) - bc(2) * bd(1)

  singlr = .true.

  do i = 4, 1, -1

    t = sqrt ( mat(i,1)**2 + mat(i,2)**2 + mat(i,3)**2 )

    if ( t <= rtol ) then
      singlr = .true.
      rad = 0.0D+00
      return
    end if

    mat(i,1:3) = mat(i,1:3) / t

    if ( i == 4 ) then
      mat(i,4) = mat(i,1) * b(1) + mat(i,2) * b(2) + mat(i,3) * b(3)
    else
      mat(i,4) = mat(i,1) * a(1) + mat(i,2) * a(2) + mat(i,3) * a(3)
      mat(i,1:4) = mat(i,1:4) - mat(4,1:4)
    end if

  end do
!
!  Use Gaussian elimination with partial pivoting to solve 3 by 3
!  system of linear equations for center of insphere.
!
  r = 1

  do i = 2, 3
    if ( abs ( mat(r,1)) < abs ( mat(i,1) ) ) then
      r = i
    end if
  end do

  if ( abs ( mat(r,1) ) <= tol ) then
    singlr = .true.
    rad = 0.0D+00
    return
  end if

  if ( r /= 1 ) then
    do j = 1, 4
      t = mat(1,j)
      mat(1,j) = mat(r,j)
      mat(r,j) = t
    end do
  end if

  do i = 2, 3
    t = mat(i,1) / mat(1,1)
    mat(i,2:4) = mat(i,2:4) - t * mat(1,2:4)
  end do

  if ( abs ( mat(2,2) ) < abs ( mat(3,2) ) ) then

    do j = 2, 4
      t = mat(2,j)
      mat(2,j) = mat(3,j)
      mat(3,j) = t
    end do

  end if

  if ( tol < abs ( mat(2,2) ) ) then

    t = mat(3,2) / mat(2,2)
    mat(3,3:4) = mat(3,3:4) - t * mat(2,3:4)

    if ( tol < abs ( mat(3,3) ) ) then
      singlr = .false.
    end if

  end if

  if ( singlr ) then
    rad = 0.0D+00
    return
  end if

  center(3) = mat(3,4) / mat(3,3)
  center(2) = (mat(2,4) - mat(2,3)*center(3))/mat(2,2)
  center(1) = (mat(1,4) - mat(1,3)*center(3) - &
    mat(1,2)*center(2))/mat(1,1)

  rad = abs ( mat(4,1) * center(1) + mat(4,2) * center(2) + &
    mat(4,3) * center(3) - mat(4,4))

  return
end
