function sangmn ( a, b, c, d, sang )

!*****************************************************************************80
!
!! SANGMN computes the minimum solid angle of a tetrahedron.
!
!  Discussion:
!
!    This routine computes the four solid angles of a tetrahedron,
!    and their minimum.
!
!    Actually, sin ( solid angle / 2 ) is computed.
!
!    Actually, the result must also be multiplied by 1.5 * sqrt ( 6 )!
!
!    In that case, the result will lie between 0 and 1, and an
!    equilateral tetrahedron will achieve a value of 1.
!
!  Modified:
!
!    16 August 2005
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
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), the vertices
!    of the tetrahedron.
!
!    Output, real ( kind = 8 ) SANG(1:4), the four solid angles at A, B,
!    C, D, respectively.  (Actually sin(angle/2)).
!
!    Output, real ( kind = 8 ) SANGMN, sin(min solid angle / 2)
!    = min(SANG(1:4)) since sum of 4 solid angles <= 2*pi.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) a(dim_num)
  real    ( kind = 8 ) ab(dim_num)
  real    ( kind = 8 ) ac(dim_num)
  real    ( kind = 8 ) ad(dim_num)
  real    ( kind = 8 ) b(dim_num)
  real    ( kind = 8 ) bc(dim_num)
  real    ( kind = 8 ) bd(dim_num)
  real    ( kind = 8 ) c(dim_num)
  real    ( kind = 8 ) cd(dim_num)
  real    ( kind = 8 ) d(dim_num)
  real    ( kind = 8 ) denom
  real    ( kind = 8 ) l1
  real    ( kind = 8 ) l2
  real    ( kind = 8 ) l3
  real    ( kind = 8 ) lab
  real    ( kind = 8 ) lac
  real    ( kind = 8 ) lad
  real    ( kind = 8 ) lbc
  real    ( kind = 8 ) lbd
  real    ( kind = 8 ) lcd
  real    ( kind = 8 ) sang(4)
  real    ( kind = 8 ) sangmn
  real    ( kind = 8 ) vol
!
!  Compute the vectors that represent the sides.
!
  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)
  ad(1:3) = d(1:3) - a(1:3)
  bc(1:3) = c(1:3) - b(1:3)
  bd(1:3) = d(1:3) - b(1:3)
  cd(1:3) = d(1:3) - c(1:3)
!
!  Compute the lengths of the sides.
!
  lab = sqrt ( sum ( ab(1:dim_num)**2 ) )
  lac = sqrt ( sum ( ac(1:dim_num)**2 ) )
  lad = sqrt ( sum ( ad(1:dim_num)**2 ) )
  lbc = sqrt ( sum ( bc(1:dim_num)**2 ) )
  lbd = sqrt ( sum ( bd(1:dim_num)**2 ) )
  lcd = sqrt ( sum ( cd(1:dim_num)**2 ) )

  vol = 2.0D+00 * abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) )

  l1 = lab + lac
  l2 = lab + lad
  l3 = lac + lad

  denom = ( l1 + lbc ) * ( l1 - lbc ) &
        * ( l2 + lbd ) * ( l2 - lbd ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    sang(1) = 0.0D+00
  else
    sang(1) = vol / sqrt ( denom )
  end if

  l1 = lab + lbc
  l2 = lab + lbd
  l3 = lbc + lbd

  denom = ( l1 + lac ) * ( l1 - lac ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    sang(2) = 0.0D+00
  else
    sang(2) = vol / sqrt ( denom )
  end if

  l1 = lac + lbc
  l2 = lac + lcd
  l3 = lbc + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lbd ) * ( l3 - lbd )

  if ( denom <= 0.0D+00 ) then
    sang(3) = 0.0D+00
  else
    sang(3) = vol / sqrt ( denom )
  end if

  l1 = lad + lbd
  l2 = lad + lcd
  l3 = lbd + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lac ) * ( l2 - lac ) &
        * ( l3 + lbc ) * ( l3 - lbc )

  if ( denom <= 0.0D+00 ) then
    sang(4) = 0.0D+00
  else
    sang(4) = vol / sqrt ( denom )
  end if

  sangmn = min ( sang(1), sang(2), sang(3), sang(4) )

  return
end
