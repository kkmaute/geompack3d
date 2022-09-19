subroutine swapes ( bndcon, i, npt, sizht, fc_num, fc_max, vcl, vm, bf, fc, ht, &
  ntetra, hdavfc, front, back, ifac, ierr )

!*****************************************************************************80
!
!! SWAPES swaps faces in a 3D triangulation.
!
!  Discussion:
!
!    This routine swaps faces, applying local transformations, in a 3D
!    triangulation based on the empty circumsphere criterion until (nearly)
!    all faces are locally optimal.  I is the index of the new vertex
!    added to the triangulation, or 0 if an initial triangulation is given.
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
!  Parameters:
!
!    Input, logical BNDCON, TRUE iff boundary faces are constrained (i.e. not
!    swapped by local transformations).
!
!    Input, integer ( kind = 4 ) I, the local index of next vertex inserted in
!    triangulation, or 0; if positive, it is assumed I is largest index so far.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the  array of boundary face records;
!    see DTRIS3.
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
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of
!    queue of interior faces for which sphere test is applied.
!
!    Output, integer ( kind = 4 ) IFAC, the index of last face for which sphere test
!    applied, or 0.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  real    ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,*)
  integer ( kind = 4 ) bfx(2)
  logical              bndcon
  integer ( kind = 4 ) c
  real    ( kind = 8 ) center(3)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  logical              degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) indx(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nbr(2,2)
  integer ( kind = 4 ) ntetra
  real    ( kind = 8 ) radsq
  real    ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  ifac = 0

  do

    do

      if ( front == 0 ) then
        return
      end if

      ind = front
      front = fc(7,ind)

      if ( fc(2,ind) /= 0 ) then
        exit
      end if

      if ( ind == fc_num ) then
        fc_num = fc_num - 1
      else
        fc(1,ind) = -hdavfc
        hdavfc = ind
      end if

    end do

    ifac = ind
    fc(7,ind) = -1
    a = fc(1,ind)
    b = fc(2,ind)
    c = fc(3,ind)
    d = fc(4,ind)
    e = fc(5,ind)
    va = vm(a)
    vb = vm(b)
    vc = vm(c)
    vd = vm(d)
    ve = vm(e)

    if ( msglvl == 4 ) then
      write ( *,600) ind,a,b,c,d,e
    end if

    call ccsph ( .true., vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), &
      vcl(1,ve), center, radsq, in )

    if ( in == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SWAPES - Fatal error!'
      write ( *, '(a)' ) '  CCSPH returned IN = 2.'
      ierr = 301
      return
    end if

    if ( 1 <= in ) then

      call baryth ( vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), &
        vcl(1,ve), alpha, degen )

      if ( degen ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SWAPES - Fatal error!'
        write ( *, '(a)' ) '  BARYTH detected a degenerate tetrahedron.'
        ierr = 301
        return
      else if ( 0.0D+00 < alpha(4) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SWAPES - Fatal error!'
        write ( *, '(a)' ) '  BARYTH detected 0 < ALPHA(4).'
        ierr = 309
        return
      end if

      kneg = 1
      kzero = 0

      do j = 1, 3
        if ( alpha(j) < 0.0D+00 ) then
          kneg = kneg + 1
        else if ( alpha(j) == 0.0D+00 ) then
          kzero = kzero + 1
        end if
      end do
!
!  Swap 2 tetrahedra for 3.
!
      if ( kneg == 1 .and. kzero == 0 ) then

        call updatf ( a, b, d, c, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, c, d, b, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( b, c, d, a, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, b, e, c, d, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, c, e, b, d, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( b, c, e, a, d, i, npt, sizht, front, back, fc, ht, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a, i6)' ) '  UPDATF returned IERR = ', ierr
          return
        end if

        call htdel ( ind, npt, sizht, fc, ht )
        call htins ( ind, a, d, e, b, c, npt, sizht, fc, ht )
        call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
          return
        end if

        call htins ( ind, b, d, e, a, c, npt, sizht, fc, ht )
        call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
          return
        end if

        call htins ( ind, c, d, e, a, b, npt, sizht, fc, ht )
        ntetra = ntetra + 1

        if ( msglvl == 4 ) then
          write ( *,610)
        end if
!
!  Swap 3 tetrahedra for 2 if possible. Relabel so edge
!  AB would be deleted. Swap if ABDE is in current triangulation.
!
      else if ( kneg == 2 .and. kzero == 0 ) then

        if ( alpha(1) < 0.0D+00 ) then
          call i4_swap ( a, c )
        else if ( alpha(2) < 0.0D+00 ) then
          call i4_swap ( b, c )
        end if

        ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

        if ( ind1 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
          ierr = 300
          return
        end if

        if ( fc(4,ind1) == e .or. fc(5,ind1) == e ) then

          call updatf ( a, c, d, b, e, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( a, c, e, b, d, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( a, d, e, b, c, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, c, d, a, e, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, c, e, a, d, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, d, e, a, c, i, npt, sizht, front, back, fc, ht, ierr )

          if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a,i6)' ) '  UPDATF returned IERR = ', ierr
            return
          end if

          call htdel ( ind, npt, sizht, fc, ht )
          call htins ( ind, c, d, e, a, b, npt, sizht, fc, ht )
          call htdel ( ind1, npt, sizht, fc, ht )

          if ( 0 <= fc(7,ind1) ) then
            fc(2,ind1) = 0
          else
            if ( ind1 == fc_num ) then
              fc_num = fc_num - 1
            else
              fc(1,ind1) = -hdavfc
              hdavfc = ind1
            end if
          end if

          ind1 = htsrc ( a, b, e, npt, sizht, fc, ht )

          if ( ind1 <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
            ierr = 300
            return
          end if

          call htdel ( ind1, npt, sizht, fc, ht )

          if ( 0 <= fc(7,ind1) ) then

            fc(2,ind1) = 0

          else

            if ( ind1 == fc_num ) then
              fc_num = fc_num - 1
            else
              fc(1,ind1) = -hdavfc
              hdavfc = ind1
            end if

          end if

          ntetra = ntetra - 1

          if ( msglvl == 4 ) then
            write ( *,620) c,d,e
          end if

        else

          if ( msglvl == 4 ) then
            write ( *,630) a,b,d,e
          end if

        end if
!
!  Coplanar faces: swap 2 tetrahedra for 2 if boundary faces
!  (and BNDCON is .FALSE.), else do pair of 2 for 2 swaps if
!  possible.  Relabel vertices so that DE intersects AB.
!  Also swap if necessary to make A < B and D < E.
!
      else if ( kneg == 1 .and. kzero == 1 ) then

        if ( alpha(1) == 0.0D+00 ) then

          call i4_swap ( a, c )

        else if ( alpha(2) == 0.0D+00 ) then

          call i4_swap ( b, c )

        end if

        if ( b < a ) then
          call i4_swap ( a, b )
        end if

        if ( e < d ) then
          call i4_swap ( d, e )
        end if

        ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

        if ( ind1 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
          ierr = 300
          return
        end if

        ind2 = htsrc ( a, b, e, npt, sizht, fc, ht )

        if ( ind2 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND2 <= 0.'
          ierr = 300
          return
        end if

        if ( fc(4,ind1) == c ) then
          f = fc(5,ind1)
        else
          f = fc(4,ind1)
        end if

        if ( fc(4,ind2) == c ) then
          g = fc(5,ind2)
        else
          g = fc(4,ind2)
        end if

        if ( f <= 0 .and. g <= 0 ) then

          if ( .not. bndcon ) then

            call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht, ierr )

            if ( ierr /= 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SWAPES - Fatal error!'
              write ( *, '(a,i6)' ) '  UPDATF returned IERR = ', ierr
              return
            end if

            call htdel(ind,npt,sizht,fc,ht)
            call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
            call htdel(ind1,npt,sizht,fc,ht)
            call htins(ind1,a,d,e,c,fc(5,ind1),npt,sizht,fc,ht)
            call htdel(ind2,npt,sizht,fc,ht)
            call htins(ind2,b,d,e,c,fc(5,ind2),npt,sizht,fc,ht)
            indx(1) = ind1
            indx(2) = ind2
            bfx(1) = -fc(5,ind1)
            bfx(2) = -fc(5,ind2)
            dd = d

            do j = 1, 2

              if ( j == 2 ) then
                dd = e
              end if

              if ( dd < a ) then
                nbr(j,1) = bf(3,bfx(j))
                nbr(j,2) = bf(2,bfx(j))
              else if ( dd < b ) then
                nbr(j,1) = bf(3,bfx(j))
                nbr(j,2) = bf(1,bfx(j))
              else
                nbr(j,1) = bf(2,bfx(j))
                nbr(j,2) = bf(1,bfx(j))
              end if

            end do

            aa = a
            k = -fc(5,nbr(1,2))

            do j = 1, 2

              if ( j == 2 ) then
                aa = b
                k = -fc(5,nbr(2,1))
              end if

              if ( aa < d ) then
                bf(1,bfx(j)) = indx(3-j)
                bf(2,bfx(j)) = nbr(2,j)
                bf(3,bfx(j)) = nbr(1,j)
              else if ( aa < e ) then
                bf(1,bfx(j)) = nbr(2,j)
                bf(2,bfx(j)) = indx(3-j)
                bf(3,bfx(j)) = nbr(1,j)
              else
                bf(1,bfx(j)) = nbr(2,j)
                bf(2,bfx(j)) = nbr(1,j)
                bf(3,bfx(j)) = indx(3-j)
              end if

              if ( bf(1,k) == indx(j) ) then
                bf(1,k) = indx(3-j)
              else if ( bf(2,k) == indx(j) ) then
                bf(2,k) = indx(3-j)
              else
                bf(3,k) = indx(3-j)
              end if

            end do

            if ( msglvl == 4 ) then
              write ( *,640) a,b,d,e
            end if

          end if

        else if ( f == g ) then

          call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,d,f,b,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,e,f,b,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,d,f,a,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,e,f,a,d,i,npt,sizht,front,back,fc,ht, ierr )

          if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a,i6)' ) '  UPDATF returned IERR = ', ierr
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)
          call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)

          ind = htsrc ( a, b, f, npt, sizht, fc, ht )

          if ( ind <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a)' ) '  HTSRC returned IND <= 0.'
            ierr = 300
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)

          if ( 0 <= fc(7,ind) ) then
            fc(2,ind) = 0
            call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )
            if ( ierr /= 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SWAPES - Fatal error!'
              write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
              return
            end if
          end if

          call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)
          call htdel(ind1,npt,sizht,fc,ht)
          j = fc(7,ind1)
          call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)
          fc(7,ind1) = j
          call htdel(ind2,npt,sizht,fc,ht)
          j = fc(7,ind2)
          call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)
          fc(7,ind2) = j

          if ( i <= 0 .and. fc(7,ind1) == -1 ) then
            fc(7,ind1) = 0
            if ( front == 0 ) then
              front = ind1
            else
              fc(7,back) = ind1
            end if
            back = ind1
          end if

          if ( i <= 0 .and. fc(7,ind2) == -1 ) then
            fc(7,ind2) = 0
            if ( front == 0 ) then
              front = ind2
            else
              fc(7,back) = ind2
            end if
            back = ind2
          end if

          if ( msglvl == 4 ) then
            write ( *,650) a,b,d,e,f
          end if

        else

          if ( msglvl == 4 ) then
            write ( *,660) a,b,d,e,f,g
          end if

        end if

      end if

    end if

  end do

  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 3-2 not poss, tetra missing:',4i7)
  640 format (4x,'swap 2-2: edge ',2i7,' repl by ',2i7)
  650 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)
  660 format (4x,'swap 4-4 not poss: a,b,d,e,f,g =',6i7)

  return
end
