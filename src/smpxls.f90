subroutine smpxls ( k, nfc, vm, fc, ns, smplx )

!*****************************************************************************80
!
!! SMPXLS constructs a list of simplices from the FC array.
!
!  Discussion:
!
!    This routine constructs a list of simplices from FC array.  Global vertex
!    indices from VM are produced.  The vertex indices for each
!    simplex are sorted in increasing order.
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
!    Input, integer ( kind = 4 ) K, the dimension of triangulation.
!
!    Input, integer ( kind = 4 ) NFC, number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) VM(1:*), indices of vertices of VCL that are triangulated.
!
!    Input, integer ( kind = 4 ) FC(1:K+4,1:NFC), array of face records; see routine DTRISK.
!
!    Output, integer ( kind = 4 ) NS, number of simplices.
!
!    Output, integer ( kind = 4 ) SMPLX(1:K+1,1:NS), contains global simplex indices; it is
!    assumed there is enough space.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) nfc

  integer ( kind = 4 ) fc(k+4,nfc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ns
  integer ( kind = 4 ) smplx(k+1,*)
  integer ( kind = 4 ) vm(*)

  kp1 = k + 1
  kp2 = k + 2
  ns = 0

  do i = 1, nfc

    if ( fc(1,i) <= 0 ) then
      cycle
    end if

    do l = kp1, kp2
      if ( fc(k,i) < fc(l,i) ) then
        ns = ns + 1
        smplx(1:k,ns) = vm(fc(1:k,i))
        smplx(kp1,ns) = vm(fc(l,i))
      end if
    end do

  end do

  do i = 1, ns
    call orderk ( kp1, smplx(1,i) )
  end do

  return
end
