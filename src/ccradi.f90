function ccradi ( a, b, c, d )

!*****************************************************************************80
!
!! CCRADI computes the circumradius of a tetrahedron.
!
!  Discussion:
!
!    This routine computes the circumradius of a tetrahedron.
!
!    Actually 1/(4*circumradius)**2 is computed.
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
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), the
!    vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) CCRADI, the value 1/(4*circumradius)**2.
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
  real    ( kind = 8 ) ccradi
  real    ( kind = 8 ) cd(3)
  real    ( kind = 8 ) d(3)
  real    ( kind = 8 ) denom
  real    ( kind = 8 ) lab
  real    ( kind = 8 ) lac
  real    ( kind = 8 ) lad
  real    ( kind = 8 ) lbc
  real    ( kind = 8 ) lbd
  real    ( kind = 8 ) lcd
  real    ( kind = 8 ) pb
  real    ( kind = 8 ) pc
  real    ( kind = 8 ) pd
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) vol

  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)
  ad(1:3) = d(1:3) - a(1:3)
  bc(1:3) = c(1:3) - b(1:3)
  bd(1:3) = d(1:3) - b(1:3)
  cd(1:3) = d(1:3) - c(1:3)

  lab = ab(1)**2 + ab(2)**2 + ab(3)**2
  lac = ac(1)**2 + ac(2)**2 + ac(3)**2
  lad = ad(1)**2 + ad(2)**2 + ad(3)**2
  lbc = bc(1)**2 + bc(2)**2 + bc(3)**2
  lbd = bd(1)**2 + bd(2)**2 + bd(3)**2
  lcd = cd(1)**2 + cd(2)**2 + cd(3)**2

  pb = sqrt ( lab * lcd )
  pc = sqrt ( lac * lbd )
  pd = sqrt ( lad * lbc )
  t1 = pb + pc
  t2 = pb - pc
  denom = (t1+pd) * (t1-pd) * (pd+t2) * (pd-t2)

  if ( denom <= 0.0D+00 ) then
    ccradi = 0.0D+00
    return
  end if

  vol = ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
      + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
      + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) )

  ccradi = vol**2 / denom

  return
end
