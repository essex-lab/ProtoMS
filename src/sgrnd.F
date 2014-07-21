      subroutine sgrnd(seed)
c
      implicit integer(a-z)
c
c Period parameters
      parameter(N     =  624)
c
      dimension mt(0:N-1)
c                     the array for the state vector
      common /block/mti,mt
      save   /block/
c
c      setting initial seeds to mt[N] using
c      the generator Line 25 of Table 1 in
c      [KNUTH 1981, The Art of Computer Programming
c         Vol. 2 (2nd Ed.), pp102]
c
      mt(0)= iand(seed,-1)
      do mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
      enddo
c
      return
      end
