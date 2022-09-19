function prime ( k )

!*****************************************************************************80
!
!! PRIME returns a prime greater than a given value K.
!
!  Discussion:
!
!    This routine returns a prime greater than or equal to K (if possible)
!    from internal array of primes.  More primes can be added if desired.
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
!    Input, integer ( kind = 4 ) K, a positive integer.
!
!    Output, integer ( kind = 4 ) PRIME, the smallest prime greater than or equal to K
!    from internal array (or largest in array).
!
  implicit none

  integer ( kind = 4 ), parameter :: nprime = 150

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) prime
  integer ( kind = 4 ), parameter, dimension ( nprime ) :: primes = (/ &
    17,31,47,61,79,97,113,127,149,163,179,193,211,227,241, &
    257,271,293,307,331,353,379,401,431,457,479,503,541,563,587, &
    613,641,673,701,727,751,773,797,821,853,877,907,929,953,977, &
    1009,1049,1087,1123,1163,1201,1237,1277,1319,1361,1399,1433, &
    1471,1511,1543,1579,1613,1657,1699,1741,1783,1831,1873,1931, &
    1973,2017,2069,2129,2203,2267,2333,2389,2441,2503,2557,2609, &
    2663,2719,2789,2851,2917,2999,3061,3137,3209,3299,3371,3449, &
    3527,3613,3697,3779,3863,3947,4049,4211,4421,4621,4813,5011, &
    5227,5413,5623,5813,6011,6211,6421,6619,6823,7013,7211,7411, &
    7621,7817,8011,8219,8419,8623,8819,9011,9221,9413,9613,9811, &
    10037,10211,10427,10613,10831,11027,11213,11411,11617,11813, &
    12011,12211,12413,12611,12821,13033,13217,13411,13613,13829, &
    14011/)
  integer ( kind = 4 ) u

  if ( k <= primes(1) ) then
    prime = primes(1)
    return
  else if ( primes(nprime) <= k ) then
    prime = primes(nprime)
    return
  end if
!
!  Use binary search to find K <= prime.
!
  l = 1
  u = nprime

  do

    m = ( l + u ) / 2

    if ( k < primes(m) ) then
      u = m - 1
    else if ( primes(m) < k ) then
      l = m + 1
    else
      prime = primes(m)
      return
    end if

    if ( u < l ) then
      exit
    end if

  end do

  prime = primes(u+1)

  return
end
