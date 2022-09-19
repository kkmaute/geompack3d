subroutine nwthou ( i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, &
  ntetra, hdavbf, hdavfc, front, back, bfi, ierr )

!*****************************************************************************80
!
!! NWTHOU creates new tetrahedra outside the current convex hull.
!
!  Discussion:
!
!    This routine creates new tetrahedra in a 3D triangulation outside the
!    convex hull by joining vertex I to visible boundary faces.
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
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:BF_MAX), the array of boundary face
!    records; see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVBF, the head pointer to available BF records.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input, integer ( kind = 4 ) FRONT, the index of front of queue (or top of stack)
!    of visible boundary faces.
!
!    Output, integer ( kind = 4 ) BACK, the index of back of queue (or bottom of stack)
!    of visible boundary faces (which become interior faces).
!
!    Output, integer ( kind = 4 ) BFI, the index of FC of a boundary face containing
!    vertex I.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfnew
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) ptr
!
!  For ABC in queue, form tetrahedron ABCI + add faces ABI, ACI, BCI.
!  PTR, NBR, IND are indices of FC; K, L, BFNEW indices of BF.
!
  ierr = 0
  bfi = 0
  ptr = front

  do

    back = ptr
    a = fc(1,ptr)
    b = fc(2,ptr)
    c = fc(3,ptr)
    k = -fc(5,ptr)
    fc(5,ptr) = i
    ntetra = ntetra + 1

    if ( msglvl == 4 ) then
      write ( *,600) a,b,c,i
    end if

    do e = 1, 3

      if ( e == 2 ) then
        call i4_swap ( a, b )
      else if ( e == 3 ) then
        call i4_swap ( a, c )
      end if

      nbr = bf(e,k)

      if ( fc(7,nbr) /= -1 ) then
        if ( fc(5,nbr) == i ) then
          cycle
        end if
      end if

      call availf ( hdavfc, nfc, fc_max, fc, ind, ierr )

      if ( ierr /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NWTHOU - Error!'
        write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
        return
      end if

      l = -fc(5,nbr)

      if ( bf(1,l) == ptr ) then
        j = 1
      else if ( bf(2,l) == ptr ) then
        j = 2
      else
        j = 3
      end if

      if ( fc(7,nbr) /= -1 ) then

        call htins ( ind, b, c, i, a, fc(j,nbr), npt, sizht, fc, ht )

      else

        if ( hdavbf /= 0 ) then
          bfnew = hdavbf
          hdavbf = -bf(1,hdavbf)
        else
          if ( bf_max <= bf_num ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'NWTHOU - Error!'
            write ( *, '(a)' ) '  BF_MAX <= BF_NUM.'
            write ( *, '(a)' ) '  Increase memory BF_MAX to proceed.'
            ierr = 12
            return
          end if
          bf_num = bf_num + 1
          bfnew = bf_num
        end if

        if ( bfi == 0 ) then
          bfi = ind
        end if

        call htins ( ind, b, c, i, a, -bfnew, npt, sizht, fc, ht )
        bf(j,l) = ind
        bf(3,bfnew) = nbr

      end if

    end do

    if ( k == bf_num ) then
      bf_num = bf_num - 1
    else
      bf(1,k) = -hdavbf
      hdavbf = k
    end if

    ptr = fc(7,ptr)

    if ( ptr == 0 ) then
      exit
    end if

  end do
!
!  Set BF(1:2,BFNEW) fields for new boundary faces.
!
  ptr = bfi
  a = fc(1,ptr)
  j = 2

  do

    b = fc(j,ptr)
    c = fc(4,ptr)

    do

      nbr = htsrc ( a, c, i, npt, sizht, fc, ht )

      if ( nbr <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NWTHOU - Error!'
        write ( *, '(a,i6)' ) '  HTSRC returned IERR = ', ierr
        ierr = 300
        return
      end if

      if ( fc(5,nbr) <= 0 ) then
        exit
      end if

      if ( fc(4,nbr) == b ) then
        d = fc(5,nbr)
      else
        d = fc(4,nbr)
      end if

      b = c
      c = d

    end do

    k = -fc(5,ptr)
    l = -fc(5,nbr)

    if ( fc(1,ptr) == a ) then
      bf(2,k) = nbr
    else
      bf(1,k) = nbr
    end if

    if ( fc(1,nbr) == a ) then
      j = 1
    else
      j = 2
    end if

    bf(3-j,l) = ptr
    a = fc(3-j,nbr)
    ptr = nbr

    if ( ptr == bfi ) then
      exit
    end if

  end do

  600 format ( '  New tetra: ',4i7)

  return
end
