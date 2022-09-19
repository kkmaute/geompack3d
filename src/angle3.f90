function angle3 ( u, v, rtolsq )

!*****************************************************************************80
!
!! ANGLE3 computes the size of a plane angle in 3D.
!
!  Discussion:
!
!    This routine computes the angle, in the range [0,PI], between two
!    3D vectors U and V.
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
!    Input, real ( kind = 8 ) U(1:3), V(1:3), the vectors.
!
!    Input, real ( kind = 8 ) RTOLSQ, the relative tolerance used to
!    detect a zero vector, based on the square of the Euclidean length.
!
!    Output, real ( kind = 8 ) ANGLE3, the angle between the two vectors,
!    in the range [0,PI].  If U or V is the zero vector, ANGLE3 = PI
!    is returned.
!
  implicit none

  real    ( kind = 8 ) angle3
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) lu
  real    ( kind = 8 ) lv
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) rtolsq
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) u(3)
  real    ( kind = 8 ) v(3)

  tol = 100.0D+00 * epsilon ( tol )

  dotp = dot_product ( u(1:3), v(1:3) )

  lu = dot_product ( u(1:3), u(1:3) )

  lv = dot_product ( v(1:3), v(1:3) )

  if ( rtolsq < lu .and. rtolsq < lv ) then
    t = dotp / sqrt ( lu * lv )
    if ( 1.0D+00 - tol < abs ( t ) ) then
      t = sign ( 1.0D+00, t )
    end if
    angle3 = acos ( t )
  else
    angle3 = pi
  end if

  return
end
