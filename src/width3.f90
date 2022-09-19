subroutine width3 ( nface, vcl, hvl, nrml, fvl, maxiw, i1, i2, wid, iwk, ierr )

!*****************************************************************************80
!
!! WIDTH3 determines the width of a convex polyhedron.
!
!  Discussion:
!
!    This routine computes the width (minimum breadth) of a convex polyhedron
!    which is stored in convex polyhedron data structure.
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
!    Input, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit outward normals
!    of faces.
!
!    Input, integer ( kind = 4 ) FVL(1:5,1:*), the face vertex list; see routine DSCPH.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be at least number of edges in polyhedron.
!
!    Output, integer ( kind = 4 ) I1, I2, the indices realizing width, either face in
!    HVL + vertex in FVL if positive, or 2 edges in FVL if negative.
!
!    Output, real ( kind = 8 ) WID, the width of convex polyhedron.
!
!    Workspace, integer IWK(1:MAXIW) - used for indices of FVL.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) nface

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) d1
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) dir1(3)
  real    ( kind = 8 ) dist
  real    ( kind = 8 ) dmax
  real    ( kind = 8 ) dotp1
  real    ( kind = 8 ) dotp2
  real    ( kind = 8 ) dtol
  integer ( kind = 4 ), parameter :: edgv = 5
  real    ( kind = 8 ) en(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) ff
  integer ( kind = 4 ) fvl(5,*)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) gg
  integer ( kind = 4 ) hvl(nface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) ld
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) ne
  real    ( kind = 8 ) nrml(3,nface)
  real    ( kind = 8 ) nrmlf(3)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) rhs
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,*)
  real    ( kind = 8 ) wid
  real    ( kind = 8 ) wtol
!
!  Determine list of distinct vertices of polyhedron, and store
!  indices of FVL for these vertices in IWK.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  nv = 0

  do f = 1,nface

    a = hvl(f)

    do

      la = fvl(loc,a)

      do i = 1, nv
        if ( la == fvl(loc,iwk(i))) go to 30
      end do

      nv = nv + 1

      if ( maxiw < nv ) then
        ierr = 6
        return
      end if

      iwk(nv) = a

30    continue

      a = fvl(succ,a)

      if ( a == hvl(f) ) then
        exit
      end if

    end do

  end do

  ne = nv + nface - 2

  if ( maxiw < ne ) then
    ierr = 6
    return
  end if
!
!  For each face, find vertex furthest from face and its distance.
!
  wid = 0.0D+00

  do f = 1,nface

    la = fvl(loc,hvl(f))
    rhs = nrml(1,f)*vcl(1,la) + nrml(2,f)*vcl(2,la) + nrml(3,f)*vcl(3,la)
    dist = 0.0D+00

    do i = 1, nv

      la = fvl(loc,iwk(i))
      d = rhs - nrml(1,f)*vcl(1,la) - nrml(2,f)*vcl(2,la) &
        - nrml(3,f)*vcl(3,la)

      if ( dist < d ) then
        dist = d
        iv = iwk(i)
      end if

    end do

    if ( dist < wid .or. f == 1 ) then
      wid = dist
      i1 = f
      i2 = iv
    end if

  end do
!
!  Determine list of edges of polyhedron, and store indices of FVL
!  for these edges in IWK.
!
  ne = 0

  do f = 1,nface

    a = hvl(f)
    lb = fvl(loc,a)

70  continue

    la = lb
    b = fvl(succ,a)
    lb = fvl(loc,b)

    if ( la < lb ) then
      ne = ne + 1
      iwk(ne) = a
    end if

    a = b

    if ( a /= hvl(f) ) then
      go to 70
    end if

  end do
!
!  For each pair of edges, find parallel planes through edges and
!  determine whether it is antipodal pair.
!
  wtol = tol * wid

  do i = 1, ne-1

    a = iwk(i)
    f = fvl(facn,a)
    ff = fvl(facn,fvl(edgv,a))
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))
    dir(1:3) = vcl(1:3,lb) - vcl(1:3,la)
    dmax = max ( abs ( dir(1) ), abs ( dir(2) ), abs ( dir(3) ) )

    do j = i+1,ne

      c = iwk(j)
      g = fvl(facn,c)
      gg = fvl(facn,fvl(edgv,c))

      if ( g == f .or. g == ff .or. gg == f .or. gg == ff ) then
        cycle
      end if

      lc = fvl(loc,c)
      ld = fvl(loc,fvl(succ,c))

      if ( lc == la .or. lc == lb .or. ld == la .or. ld == lb)  then
        cycle
      end if

      dir1(1:3) = vcl(1:3,ld) - vcl(1:3,lc)
      dtol = tol * max ( dmax, abs ( dir1(1) ), abs ( dir1(2) ), &
        abs ( dir1(3) ) )
      nrmlf(1) = dir(2)*dir1(3) - dir(3)*dir1(2)
      nrmlf(2) = dir(3)*dir1(1) - dir(1)*dir1(3)
      nrmlf(3) = dir(1)*dir1(2) - dir(2)*dir1(1)
      k = 1
      if ( abs(nrmlf(1)) < abs(nrmlf(2))) k = 2
      if ( abs(nrmlf(k)) < abs(nrmlf(3))) k = 3

      if ( abs(nrmlf(k)) <= dtol ) then
        cycle
      end if

      leng = sqrt(nrmlf(1)**2 + nrmlf(2)**2 + nrmlf(3)**2)
      d = nrmlf(1)*vcl(1,la) + nrmlf(2)*vcl(2,la) + &
           nrmlf(3)*vcl(3,la)
      d1 = nrmlf(1)*vcl(1,lc) + nrmlf(2)*vcl(2,lc) + &
           nrmlf(3)*vcl(3,lc)
      dist = abs(d - d1)/leng

      if ( dist <= wtol .or. wid - wtol <= dist ) then
        cycle
      end if

      nrmlf(1:3) = nrmlf(1:3) / leng

      en(1) = nrml(2,f)*dir(3) - nrml(3,f)*dir(2)
      en(2) = nrml(3,f)*dir(1) - nrml(1,f)*dir(3)
      en(3) = nrml(1,f)*dir(2) - nrml(2,f)*dir(1)
      dotp1 = (en(1)*nrmlf(1) + en(2)*nrmlf(2) + en(3)*nrmlf(3)) &
        /sqrt(en(1)**2 + en(2)**2 + en(3)**2)
      if ( abs(dotp1) <= tol) dotp1 = 0.0D+00
      en(1) = nrml(2,ff)*dir(3) - nrml(3,ff)*dir(2)
      en(2) = nrml(3,ff)*dir(1) - nrml(1,ff)*dir(3)
      en(3) = nrml(1,ff)*dir(2) - nrml(2,ff)*dir(1)
      dotp2 = -(en(1)*nrmlf(1) + en(2)*nrmlf(2) + en(3)*nrmlf(3)) &
        /sqrt(en(1)**2 + en(2)**2 + en(3)**2)

      if ( abs(dotp2) <= tol) dotp2 = 0.0D+00

      if ( dotp1*dotp2 < 0.0D+00 ) then
        cycle
      end if

      en(1) = nrml(2,g)*dir1(3) - nrml(3,g)*dir1(2)
      en(2) = nrml(3,g)*dir1(1) - nrml(1,g)*dir1(3)
      en(3) = nrml(1,g)*dir1(2) - nrml(2,g)*dir1(1)
      dotp1 = (en(1)*nrmlf(1) + en(2)*nrmlf(2) + en(3)*nrmlf(3)) &
        /sqrt(en(1)**2 + en(2)**2 + en(3)**2)
      if ( abs(dotp1) <= tol) dotp1 = 0.0D+00
      en(1) = nrml(2,gg)*dir1(3) - nrml(3,gg)*dir1(2)
      en(2) = nrml(3,gg)*dir1(1) - nrml(1,gg)*dir1(3)
      en(3) = nrml(1,gg)*dir1(2) - nrml(2,gg)*dir1(1)
      dotp2 = -(en(1)*nrmlf(1) + en(2)*nrmlf(2) + en(3)*nrmlf(3)) &
        /sqrt(en(1)**2 + en(2)**2 + en(3)**2)

      if ( abs(dotp2) <= tol) dotp2 = 0.0D+00

      if ( 0.0D+00 <= dotp1 * dotp2 ) then
        wid = dist
        i1 = -a
        i2 = -c
      end if

    end do

  end do

  return
end
