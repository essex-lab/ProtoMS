      subroutine calculateCLJPairs
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'flags.inc'
c#################################################
c
c     This subroutine runs through all of the
c     pairs of parameter types used in the simulation
c     and calculates the Chg(i)*Chg(j), Aij and
c     Bij for these parameters...
c
c     (C) Christopher Woods
c
c#################################################

      integer i,j,ncljpairs
      double precision epspair,sigpair
      double precision charge,sigma,epslon

      write(printstring,10) NCLJTypes
10    format("There are ",i4," CLJ parameter types used in this simulation.")
      call printLine(INFO,printstring)

      write(printstring,*) "Pre-calculating all the charge and LJ pairs..."
      call printLine(INFO,printstring)

      ncljpairs = 0
      do i=1,NCLJTypes
        do j=i,NCLJTypes
          ncljpairs = ncljpairs + 1
          if (ncljpairs.gt.MAXCLJPAIRS) then
            write(printstring,11) MAXCLJPAIRS
11          format("Exceeded maximum number of CLJ atom pairs(",i4,")")
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXCLJPAIRS and recompile")
            return
          endif
          
c         save the pair index...
          CLJPair(i,j) = ncljpairs
          CLJPair(j,i) = ncljpairs
          
c         calculate Chg(i) * Chg(j) / 4_pi_eps0
c         (units of kcal mol-1) when distances in A 
          ChgIJ(ncljpairs) = charge(i)*charge(j) * FAC4PIEPS0
          
          if (lorentzberthelot) then
c           calculate Aij and Bij using Lorentz-Berthelot combining rules...
c                   Aij = 4 * sqrt(eps(i)eps(j)) * (0.5*(sig(i)+sig(j)))**12
c                   Bij = 4 * sqrt(eps(i)eps(j)) * (0.5*(sig(i)+sig(j)))**6
c           Again, give units of kcal mol-1 when distances in A
            epspair = FOUR * dsqrt(epslon(i)*epslon(j))
            sigpair = (HALF*(sigma(i)+sigma(j)))**6
            AIJ(ncljpairs) = epspair * sigpair**2
            BIJ(ncljpairs) = epspair * sigpair
          else
c           calculate Aij and Bij using OPLS combining rules...
c                   Aij = 4 * sqrt(eps(i)eps(j)) * (sqrt(sig(i)sig(j)))**12
c                   Bij = 4 * sqrt(eps(i)eps(j)) * (sqrt(sig(i)sig(j)))**6
c           Again, give units of kcal mol-1 when distances in A
            epspair = FOUR * dsqrt(epslon(i)*epslon(j))
            sigpair = (sigma(i)*sigma(j))**3
            AIJ(ncljpairs) = epspair * sigpair**2
            BIJ(ncljpairs) = epspair * sigpair
          endif
        enddo
      enddo

      return 
      
      end
