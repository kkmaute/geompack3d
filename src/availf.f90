subroutine availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

!*****************************************************************************80
!
!! AVAILF returns the index of the next available record in the FC array.
!
!  Discussion:
!
!    This routine returns the index of the next available record in the
!    FC array, either HDAVFC or FC_NUM+1.
!
!  Modified:
!
!    07 September 2005
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
!    Input/output, integer ( kind = 4 ) HDAVFC, head pointer of available records in FC.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, current number of records used in FC.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum number of records available in FC.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Output, integer ( kind = 4 ) IND, the index of available record (if FC not full).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) fc_max

  ierr = 0

  if ( hdavfc /= 0 ) then
    ind = hdavfc
    hdavfc = -fc(1,hdavfc)
  else if ( fc_max <= fc_num ) then
    ierr = 11
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVAILF - Fatal error!'
    write ( *, '(a)' ) '  Memory requirements for array FC exceed the'
    write ( *, '(a,i12)' ) '  current limit of FC_MAX = ', fc_max
  else
    fc_num = fc_num + 1
    ind = fc_num
  end if

  return
end
