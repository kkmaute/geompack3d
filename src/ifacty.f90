subroutine ifacty ( ind, npt, sizht, vcl, vm, fc, ht, type, a, b, c, ierr )

!*****************************************************************************80
!
!! IFACTY determines the type of an interior face in a 3D triangulation.
!
!  Discussion:
!
!    This routine determines the type of an interior face of a 3D triangulation.
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
!    Input, integer ( kind = 4 ) IND, index of FC array, assumed to be an interior face.
!
!    Input, integer ( kind = 4 ) NPT, size of VM array.
!
!    Input, integer ( kind = 4 ) SIZHT, size of HT array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), vertex mapping array, from local to
!    global indices.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:SIZHT-1), hash table using direct chaining.
!
!    Output, character ( len = 3 ) TYPE, 'T23', 'T32', 'T22', 'T44', 'N32',
!    'N44', 'N40', 'N30',or 'N20'.
!
!    Output, integer ( kind = 4 ) A, B, C, local indices of interior face; for T32, N32,
!    T22, T44, or N44 face, AB is edge that would get swapped out and C
!    is third vertex; for T23 face, A < B < C; for N40 face,
!    A is interior vertex; for N30 face, A is inside a face
!    with vertex B; for N20 face, A is on an edge.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer   ( kind = 4 ) npt
  integer   ( kind = 4 ) sizht

  integer   ( kind = 4 ) a
  real      ( kind = 8 ) alpha(4)
  integer   ( kind = 4 ) b
  integer   ( kind = 4 ) c
  integer   ( kind = 4 ) d
  logical                degen
  integer   ( kind = 4 ) e
  integer   ( kind = 4 ) f
  integer   ( kind = 4 ) fc(7,*)
  integer   ( kind = 4 ) ht(0:sizht-1)
  integer   ( kind = 4 ) htsrc
  integer   ( kind = 4 ) ierr
  integer   ( kind = 4 ) ind
  integer   ( kind = 4 ) ind1
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) kneg
  integer   ( kind = 4 ) kzero
  character ( len = 3 )  type
  real      ( kind = 8 ) vcl(3,*)
  integer   ( kind = 4 ) vm(npt)

  ierr = 0
  a = fc(1,ind)
  b = fc(2,ind)
  c = fc(3,ind)
  d = fc(4,ind)
  e = fc(5,ind)

  call baryth(vcl(1,vm(a)),vcl(1,vm(b)),vcl(1,vm(c)),vcl(1,vm(d)), &
    vcl(1,vm(e)),alpha,degen)

  if ( degen ) then
    ierr = 301
    return
  else if ( 0.0D+00 < alpha(4) ) then
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

  type = 'xxx'

  if ( kneg == 1 .and. kzero == 0 ) then

    type = 'T23'

  else if ( kneg == 2 .and. kzero == 0 ) then

    if ( alpha(1) < 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) < 0.0D+00 ) then
      call i4_swap ( b, c )
    end if

    ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

    if ( ind1 <= 0 ) then
      ierr = 300
      return
    else if ( fc(4,ind1) == e .or. fc(5,ind1) == e ) then
      type = 'T32'
    else
      type = 'N32'
    end if

  else if ( kneg == 3 .and. kzero == 0 ) then

    type = 'N40'

    if ( 0.0D+00 < alpha(2) ) then
      call i4_swap ( a, b )
    else if ( 0.0D+00 < alpha(3) ) then
      call i4_swap ( a, c )
    end if

  else if ( kneg == 1 .and. kzero == 1 ) then

    if ( alpha(1) == 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) == 0.0D+00 ) then
      call i4_swap ( b, c )
    end if

    ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

    if ( ind1 <= 0 ) then

      ierr = 300
      return

    else if ( fc(5,ind1) <= 0 ) then

      type = 'T22'

    else

      if ( fc(4,ind1) == c ) then
        f = fc(5,ind1)
      else
        f = fc(4,ind1)
      end if

      ind1 = htsrc(a,b,e,npt,sizht,fc,ht)

      if ( ind1 <= 0 ) then
        ierr = 300
        return
      else if ( fc(4,ind1) == f .or. fc(5,ind1) == f ) then
        type = 'T44'
      else
        type = 'N44'
      end if

    end if

  else if ( kneg == 2 .and. kzero == 1 ) then

    type = 'N30'

    if ( alpha(1) == 0.0D+00 ) then
      call i4_swap ( a, c )
    else if ( alpha(2) == 0.0D+00 ) then
      call i4_swap ( b, c )
      alpha(2) = alpha(3)
    end if

    if ( 0.0D+00 < alpha(2) ) then
      call i4_swap ( a, b )
    end if

  else if ( kneg == 1 .and. kzero == 2 ) then

    type = 'N20'

    if ( 0.0D+00 < alpha(2) ) then
      call i4_swap ( a, b )
    else if ( 0.0D+00 < alpha(3) ) then
      call i4_swap ( a, c )
     end if

  end if

  return
end
