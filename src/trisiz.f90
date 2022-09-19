subroutine trisiz ( ntrid, npolg, hvl, pvl, area, psi, h, indp, loch )

!*****************************************************************************80
!
!! TRISIZ smooths PSI and computes triangle sizes.
!
!  Discussion:
!
!    This routine smooths PSI (mean mesh distribution function) values using
!    a heap so that they differ by a factor of at most 4 in adjacent
!    polygons and then computes triangle sizes for each polygon.
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
!    Input, integer ( kind = 4 ) NTRID, the desired number of triangles in mesh.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of polygons or positions used in
!    HVL array.
!
!    Input, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), the polygon vertex list.
!
!    Input, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons
!    in decomposition.
!
!    Input/output, real ( kind = 8 ) PSI(1:NPOLG), the mean mdf values in
!    the convex polygons.  On output, values are smoothed.
!
!    Output, real ( kind = 8 ) H(1:NPOLG), the triangle size for convex
!    polygons.
!
!    Workspace, integer INDP(1:NPOLG), the indices of polygon or PSI which
!    are maintained in heap according to PSI values.
!
!    Workspace, integer LOCH(1:NPOLG), the location of polygon indices in heap.
!
  implicit none

  integer ( kind = 4 ) npolg

  real    ( kind = 8 ) area(npolg)
  integer ( kind = 4 ), parameter :: edgv = 4
  real    ( kind = 8 ) factor
  real    ( kind = 8 ) h(npolg)
  integer ( kind = 4 ) hvl(npolg)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(npolg)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) loch(npolg)
  integer ( kind = 4 ) ntrid
  integer ( kind = 4 ), parameter :: polg = 2
  real    ( kind = 8 ) psi(npolg)
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ) r
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sum2

  factor = 0.25D+00

  do i = 1, npolg
    indp(i) = i
    loch(i) = i
  end do

  k = int(npolg/2)

  do l = k, 1, -1
    call sfdwmf ( l, npolg, psi, indp, loch )
  end do

  do r = npolg, 2, -1

    j = indp(1)
    indp(1) = indp(r)
    loch(indp(1)) = 1
    call sfdwmf ( 1, r-1, psi, indp, loch )
    i = hvl(j)

    do

      k = pvl(edgv,i)

      if ( 0 < k ) then
        k = pvl(polg,k)
        if ( psi(k) < psi(j) * factor ) then
          psi(k) = psi(j) * factor
          call sfupmf ( loch(k), psi, indp, loch )
        end if
      end if

      i = pvl(succ,i)

      if ( i == hvl(j) ) then
        exit
      end if

    end do

  end do

  sum2 = dot_product ( psi(1:npolg), area(1:npolg) )

  factor = 2.0D+00 / real ( ntrid, kind = 8 )

  psi(1:npolg) = psi(1:npolg) / sum2

  h(1:npolg) = sqrt ( factor / psi(1:npolg) )

  return
end
