      double precision function grnd()
c
      implicit integer(a-z)
c
c Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
c                                    constant vector a
      parameter(UMASK = -2147483647)
c                                    most significant w-r bits
      parameter(LMASK =  2147483647)
c                                    least significant r bits
c Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
c
      dimension mt(0:N-1)
c                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
c                     mti==N+1 means mt[N] is not initialized
c
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
c                        mag01(x) = x * MATA for x=0,1
c
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
c
      if(mti.ge.N) then
c                       generate N words at one time
        if(mti.eq.N+1) then
c                            if sgrnd() has not been called,
          call sgrnd(4357)
c                              a default initial seed is used
        endif
c
        do kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        do kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
c
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
c
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
c
      return
      end
