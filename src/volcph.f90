function volcph ( nface, vcl, hvl, fvl )

!*****************************************************************************80
!
!! VOLCPH computes the volume of a convex polyhedron.
!
!  Discussion:
!
!    This routine computes the volume of a convex polyhedron which is stored in
!    a convex polyhedron data structure.
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
!    Input, integer ( kind = 4 ) NFACE, the number of faces in convex polyhedron.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NFACE), the head vertex list.
!
!    Input, integer ( kind = 4 ) FVL(1:5,1:*), the face vertex list; see routine DSCPH.
!
!    Output, real ( kind = 8 ) VOLCPH, the volume of convex polyhedron.
!
  implicit none

  integer ( kind = 4 ) nface

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) cntrf(3)
  integer ( kind = 4 ), parameter :: edgv = 5
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(5,*)
  integer ( kind = 4 ) hvl(nface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) vcl(3,*)
  real    ( kind = 8 ) vol
  real    ( kind = 8 ) volcph
  real    ( kind = 8 ) volth
!
!  Compute point CNTR in polyhedron by taking weighted average of
!  vertex coordinates; weight depends on number of occurrences of vertex.
!
  n = 0
  cntr(1:3) = 0.0D+00

  do i = 1, nface

    a = hvl(i)

    do

      la = fvl(loc,a)
      n = n + 1
      cntr(1:3) = cntr(1:3) + vcl(1:3,la)
      a = fvl(succ,a)

      if ( a == hvl(i) ) then
        exit
      end if

    end do

  end do

  cntr(1:3) = cntr(1:3) / real ( n, kind = 8 )
!
!  Use CNTR to form tetrahedra with each face, and sum up volume
!  of tetrahedra.
!
  vol = 0.0D+00

  do i = 1, nface

    n = 0
    cntrf(1:3) = 0.0D+00
    a = hvl(i)

30  continue

    la = fvl(loc,a)
    n = n + 1
    cntrf(1:3) = cntrf(1:3) + vcl(1:3,la)
    a = fvl(succ,a)

    if ( a /= hvl(i) ) then
      go to 30
    end if

    cntrf(1:3) = cntrf(1:3) / real ( n, kind = 8 )
    lb = fvl(loc,a)

40  continue

    la = lb
    b = fvl(succ,a)
    lb = fvl(loc,b)
    vol = vol + volth ( cntr, cntrf, vcl(1,la), vcl(1,lb) )
    a = b

    if ( a /= hvl(i) ) then
      go to 40
    end if

  end do

  volcph = vol/6.0D+00

  return
end
