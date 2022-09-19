subroutine imptr3 ( bndcon, postlt, crit, npt, sizht, fc_max, vcl, vm, nfc, &
  ntetra, bf, fc, ht, nface, ierr )

!*****************************************************************************80
!
!! IMPTR3 improves a 3D triangulation.
!
!  Discussion:
!
!    This routine improves a given 3D triangulation by applying local
!    transformations based on some local criterion.
!
!    BF, FC, HT should be as output by DTRIS3 or DTRIW3.
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
!    Input, logical BNDCON, TRUE iff boundary faces are constrained (i.e. not
!    swapped by local transformations).
!
!    Input, logical POSTLT, TRUE iff further local transformations applied by
!    postprocessing routine IMPTRF.
!
!    Input, integer ( kind = 4 ) CRIT, the criterion code: 1 for (local max-min) solid angle
!    criterion, 2 for radius ratio criterion, 3 for mean ratio
!    criterion, any other value for empty circumsphere criterion.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices (points).
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT; a good choice is a
!    prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the array of boundary face records;
!    see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation; NFACE <= NFC.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,*)
  logical              bndcon
  integer ( kind = 4 ) crit
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  logical              postlt
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vm(npt)
!
!  Create initial queue of interior faces.
!
  ierr = 0
  hdavbf = fc(7,1)
  hdavfc = fc(7,2)
  fc(7,1) = -1
  fc(7,2) = -1
  front = 0

  do i = 1, nfc
    if ( 0 < fc(1,i) .and. 0 < fc(5,i) ) then
      if ( front == 0 ) then
        front = i
      else
        fc(7,back) = i
      end if
      back = i
    end if
  end do

  if ( front /= 0 ) then
    fc(7,back) = 0
  end if

  if ( msglvl == 4 ) then
    write ( *,600) crit
  end if

  if ( 1 <= crit .and. crit <= 3 ) then
    call swapmu(bndcon,crit,npt,sizht,nfc,fc_max,vcl,vm,bf,fc,ht, &
      ntetra,hdavfc,front,back,i, ierr )
  else
    call swapes ( bndcon, 0, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
      ntetra, hdavfc, front, back, i, ierr )
  end if

  if ( ierr /= 0 ) then
    return
  end if

  if ( 1 <= crit .and. crit <= 3 .and. postlt ) then
    call imptrf(bndcon,crit,npt,sizht,fc_max,vcl,vm,nfc,ntetra, &
      hdavfc,bf,fc,ht,ierr)
  else if ( postlt ) then
    call imptrd(bndcon,npt,sizht,fc_max,vcl,vm,nfc,ntetra,hdavfc, &
      bf,fc,ht,ierr)
  end if

  if ( ierr /= 0 ) then
    return
  end if

  nface = nfc
  ptr = hdavfc

  do while ( ptr /= 0 )
    nface = nface - 1
    ptr = -fc(1,ptr)
  end do

  fc(7,1) = hdavbf
  fc(7,2) = hdavfc

  600 format (/1x,'imptr3: criterion =',i3)

  return
end
