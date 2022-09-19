subroutine tetlst ( nfc, vm, fc, nt, tetra )

!*****************************************************************************80
!
!! TETLST constructs a list of tetrahedra from the FC array.
!
!  Discussion:
!
!    This routine constructs a list of tetrahedra from the FC array.
!
!    Global vertex indices from VM are produced.  The vertex indices for each
!    tetrahedron are sorted in increasing order.
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
!    Input, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL that are
!    triangulated.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:NFC), array of face records; see routine DTRIS3.
!
!    Output, integer ( kind = 4 ) NT, the number of tetrahedra.
!
!    Output, integer ( kind = 4 ) TETRA(1:4,1:NT), contains global tetrahedron indices; it
!    is assumed there is enough space.
!
  implicit none

  integer ( kind = 4 ) nfc

  integer ( kind = 4 ) a
  integer ( kind = 4 ) fc(7,nfc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) t(4)
  integer ( kind = 4 ) tetra(4,*)
  integer ( kind = 4 ) vm(*)

  nt = 0

  do i = 1, nfc

    if ( fc(1,i) <= 0 )  then
      cycle
    end if

    do k = 4, 5
      if ( fc(3,i) < fc(k,i) ) then
        nt = nt + 1
        tetra(1,nt) = fc(1,i)
        tetra(2,nt) = fc(2,i)
        tetra(3,nt) = fc(3,i)
        tetra(4,nt) = fc(k,i)
      end if
    end do

  end do

  do k = 1, nt

    t(1:4) = vm(tetra(1:4,k))

    do i = 1, 3
      l = i
      do j = i+1, 4
        if ( t(j) < t(l) ) then
          l = j
        end if
      end do
      a = t(i)
      t(i) = t(l)
      t(l) = a
    end do

    tetra(1:4,k) = t(1:4)

  end do

  return
end
