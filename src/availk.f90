subroutine availk ( k, hdavfc, nfc, fc_max, fc, pos, ierr )

!*****************************************************************************80
!
!! AVAILK returns the position of the next available record in the FC array.
!
!  Discussion:
!
!    This routine returns the position of the next available record in
!    the FC array, either HDAVFC or NFC+1.
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
!    Input, integer ( kind = 4 ) K, the number of vertices in a face.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, head pointer of available records in FC.
!
!    Input/output, integer ( kind = 4 ) NFC, current number of records used in FC.
!
!    Input, integer ( kind = 4 ) FC_MAX, maximum number of records available in FC.
!
!    Input, integer ( kind = 4 ) FC(1:K+4,1:*), array of face records; see routine DTRISK.
!
!    Output, integer ( kind = 4 ) POS, position of available record (if FC not full).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) pos

  ierr = 0

  if ( hdavfc /= 0 ) then
    pos = hdavfc
    hdavfc = -fc(1,hdavfc)
  else if ( fc_max <= nfc ) then
    ierr = 22
  else
    nfc = nfc + 1
    pos = nfc
  end if

  return
end
