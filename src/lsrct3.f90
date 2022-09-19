subroutine lsrct3 ( pt, n, p, nfc, vcl, vm, fc, ht, ifac, ivrt, ierr )

!*****************************************************************************80
!
!! LSRCT3 searches a 3D triangulation for the tetrahedron containing a point.
!
!  Discussion:
!
!    This routine performs a linear search through a 3D triangulation to find
!    a tetrahedron containing the point PT.
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
!    Input, real ( kind = 8 ) PT(1:3), a 3D point.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices or size of VM.
!
!    Input, integer ( kind = 4 ) P, the size of the hash table.
!
!    Input, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:N), the vertex mapping list (maps from local indices
!    used in FC to indices of VCL).
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) IFAC, the index of FC indicating tetrahedron or face
!    containing PT or 0 if PT outside convex hull.
!
!    Output, integer ( kind = 4 ) IVRT; 4 or 5 to indicate that FC(IVRT,IFAC) is 4th
!    vertex of tetrahedron containing PT in its interior; 6 if PT lies
!    in interior of face FC(*,IFAC); 1, 2, or 3 if PT lies on
!    interior of edge of face from vertices FC(IVRT,IFAC) to
!    FC(IVRT mod 3 + 1,IFAC); -1, -2, or -3 if PT is (nearly)
!    vertex FC(-IVRT,IFAC); 0 if PT lies outside convex hull.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  logical              degen
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nfc
  real    ( kind = 8 ) pt(3)
  real    ( kind = 8 ) t(4)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) vm(n)
  logical              zero(4)

  ierr = 0

  do f = 1, nfc

    a = fc(1,f)

    if ( a <= 0 ) then
      cycle
    end if

    b = fc(2,f)
    c = fc(3,f)
    ifac = f

    do i = 4, 5

      d = fc(i,f)

      if ( d <= c ) then
        cycle
      end if

      va = vm(a)
      vb = vm(b)
      vc = vm(c)
      vd = vm(d)

      call baryth ( vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), pt, t, &
        degen )

      if ( degen ) then
        ierr = 301
        return
      end if

      if ( 0.0D+00 < t(1) .and. &
           0.0D+00 < t(2) .and. &
           0.0D+00 < t(3) .and. &
           0.0D+00 < t(4) ) then
        ivrt = i
        return
      else if ( t(1) < 0.0D+00 .or. t(2) < 0.0D+00 .or. &
        t(3) < 0.0D+00 .or. t(4) < 0.0D+00 ) then
        cycle
      end if
!
!  All 0.0 <= T(J) and at least one T(J) == 0.0D+00
!
      k = 0
      do j = 1,4
        zero(j) = ( t(j) == 0.0D+00 )
        if ( zero(j) ) then
          k = k + 1
        end if
      end do

      if ( k == 1 ) then

        ivrt = 6

        if ( zero(1) ) then

          ifac = htsrc(b,c,d,n,p,fc,ht)

          if ( ifac <= 0) then
            ierr = 300
            return
          end if

        else if ( zero(2) ) then

          ifac = htsrc(a,c,d,n,p,fc,ht)

          if ( ifac <= 0) then
            ierr = 300
            return
          end if

        else if ( zero(3) ) then

          ifac = htsrc(a,b,d,n,p,fc,ht)

          if ( ifac <= 0 ) then
            ierr = 300
            return
          end if

        end if

      else if ( k == 2 ) then

        if ( zero(4) ) then

          if ( zero(3) ) then
            ivrt = 1
          else if ( zero(1) ) then
            ivrt = 2
          else
            ivrt = 3
          end if

        else

          if ( zero(3) ) then

            ifac = htsrc(a,b,d,n,p,fc,ht)

            if ( ifac <= 0) then
              ierr = 300
              return
            end if

            if ( zero(2) ) then
              aa = a
            else
              aa = b
            end if

          else

            ifac = htsrc(a,c,d,n,p,fc,ht)

            if ( ifac <= 0 ) then
              ierr = 300
              return
            end if

            aa = c
          end if

          bb = d

          if ( bb < aa ) then
            call i4_swap ( aa, bb )
          end if

          if ( fc(1,ifac) == aa ) then
            if ( fc(2,ifac) == bb ) then
              ivrt = 1
            else
              ivrt = 3
            end if
          else
            ivrt = 2
          end if

        end if

      else
!
!  K == 3
!
        if ( .not. zero(1) ) then

          ivrt = -1

        else if ( .not. zero(2) ) then

          ivrt = -2

        else if ( .not. zero(3) ) then

          ivrt = -3

        else

          ifac = htsrc(a,b,d,n,p,fc,ht)

          if ( ifac <= 0 ) then
            ierr = 300
            return
          end if

          if ( fc(1,ifac) == d ) then
            ivrt = -1
          else if ( fc(2,ifac) == d ) then
            ivrt = -2
          else
            ivrt = -3
          end if

        end if

      end if

      return

    end do

  end do

  ifac = 0
  ivrt = 0

  return
end
