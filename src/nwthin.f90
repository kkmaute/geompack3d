subroutine nwthin ( i, ifac, ivrt, npt, sizht, nfc, fc_max, fc, ht, ntetra, &
  hdavfc, front, back, ierr )

!*****************************************************************************80
!
!! NWTHIN creates new tetrahedra after the insertion of an interior vertex.
!
!  Discussion:
!
!    This routine creates new tetrahedra in a 3D triangulation from the
!    insertion of vertex I in interior of tetrahedron.
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
!    Input, integer ( kind = 4 ) I, the (local) index of next vertex inserted in
!    triangulation.
!
!    Input, integer ( kind = 4 ) IFAC, the face of tetrahedron containing vertex I
!    is FC(*,IFAC).
!
!    Input, integer ( kind = 4 ) IVRT, the 4 or 5 where 4th vertex of tetrahedron
!    is FC(IVRT,IFAC).
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of queue
!    of interior faces ABC such that ABCI is a new tetrahedron.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) indx(6)
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra

  ierr = 0
  front = 0
  back = 0
  ntetra = ntetra + 3
  a = fc(1,ifac)
  b = fc(2,ifac)
  c = fc(3,ifac)
  d = fc(ivrt,ifac)

  do j = 1, 4

    if ( j == 1 ) then

      aa = a
      bb = b
      cc = c
      dd = d
      ind = ifac

    else

      if ( j == 2 ) then
        cc = d
        dd = c
      else if ( j == 3 ) then
        bb = c
        dd = b
      else
        aa = b
        dd = a
      end if

      ind = htsrc(aa,bb,cc,npt,sizht,fc,ht)

      if ( ind <= 0 ) then
        ierr = 300
        return
      end if

    end if

    if ( fc(4,ind) == dd ) then
      fc(4,ind) = i
    else
      fc(5,ind) = i
    end if

    if ( 0 < fc(5,ind) ) then
      if ( front == 0 ) then
        front = ind
      else
        fc(7,back) = ind
      end if
      back = ind
    end if

    if ( msglvl == 4 ) then
      write ( *,600) aa,bb,cc,i
    end if

  end do

  if ( front /= 0 ) then
    fc(7,back) = 0
  end if

  do j = 1, 6
    call availf ( hdavfc, nfc, fc_max, fc, indx(j), ierr )
    if ( ierr /= 0 ) then
      return
    end if
  end do

  call htins(indx(1),a,b,i,c,d,npt,sizht,fc,ht)
  call htins(indx(2),a,c,i,b,d,npt,sizht,fc,ht)
  call htins(indx(3),a,d,i,b,c,npt,sizht,fc,ht)
  call htins(indx(4),b,c,i,a,d,npt,sizht,fc,ht)
  call htins(indx(5),b,d,i,a,c,npt,sizht,fc,ht)
  call htins(indx(6),c,d,i,a,b,npt,sizht,fc,ht)

  600 format ( '  New tetra: ',4i7)

  return
end
