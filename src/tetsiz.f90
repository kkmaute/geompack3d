subroutine tetsiz ( ntetd, npolh, facep, hfl, pfl, vol, psi, h, indp, loch )

!*****************************************************************************80
!
!! TETSIZ smooths PSI values and computes tetrahedron sizes.
!
!  Discussion:
!
!    This routine smooths PSI (mean mesh distribution function) values using
!    a heap so that they differ by a factor of at most 8 in adjacent
!    polyhedra, then computes tetrahedron sizes for each polyhedron.
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
!    Input, integer ( kind = 4 ) NTETD, the desired number of tetrahedra in mesh.
!
!    Input, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions used in
!    HFL array.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices in PFL
!    for each polyhedron.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:*), the list of signed face indices for
!    each polyhedron.
!
!    Input, real ( kind = 8 ) VOL(1:NPOLH), the volume of convex polyhedra
!    in decomposition.
!
!    Input/output, real ( kind = 8 ) PSI(1:NPOLH), the mean mdf values in
!    the convex polyhedra.  On output, the values have been smoothed.
!
!    Output, real ( kind = 8 ) H(1:NPOLH), the tetrahedron size for convex
!    polyhedra.
!
!    Workspace, integer INDP(1:NPOLH), the indices of polyhedron or PSI
!    which are maintained in heap according to PSI values.
!
!    Workspace, integer LOCH(1:NPOLH), the location of polyhedron indices
!    in heap.
!
  implicit none

  integer ( kind = 4 ) npolh

  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  real    ( kind = 8 ) factor
  real    ( kind = 8 ) h(npolh)
  integer ( kind = 4 ) hfl(npolh)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(npolh)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loch(npolh)
  integer ( kind = 4 ) ntetd
  integer ( kind = 4 ) pfl(2,*)
  real    ( kind = 8 ) power
  real    ( kind = 8 ) psi(npolh)
  integer ( kind = 4 ) r
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) vol(npolh)

  factor = 0.125D+00

  do i = 1, npolh
    indp(i) = i
    loch(i) = i
  end do

  k = int ( npolh / 2 )

  do l = k, 1, -1
    call sfdwmf ( l, npolh, psi, indp, loch )
  end do

  do r = npolh, 2, -1

    j = indp(1)
    indp(1) = indp(r)
    loch(indp(1)) = 1
    call sfdwmf ( 1, r-1, psi, indp, loch )
    i = hfl(j)

    do

      f = abs ( pfl(1,i) )

      if ( abs ( facep(2,f) ) == j ) then
        k = abs ( facep(3,f) )
      else
        k = abs ( facep(2,f) )
      end if

      if ( 0 < k ) then

        if ( psi(k) < psi(j) * factor ) then
          psi(k) = psi(j) * factor
          call sfupmf ( loch(k), psi, indp, loch )
        end if

      end if

      i = pfl(2,i)

      if ( i == hfl(j) ) then
        exit
      end if

    end do

  end do

  sum2 = dot_product ( psi(1:npolh), vol(1:npolh) )

  factor = 6.0D+00 / real ( ntetd, kind = 8)
  power = 1.0D+00 / 3.0D+00

  do i = 1, npolh
    psi(i) = psi(i) / sum2
    h(i) = ( factor / psi(i) )**power
  end do

  return
end
