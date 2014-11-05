      logical function solThetaMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'stacks.inc'
      include 'solvents.inc'
      include 'metropolis.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'moves.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'simparams.inc'
c###############################################################
c
c     This function performs everything necessary for a 
c     successful theta-solute theta move. 
c
c     (C) Samuel Genheden, from thetaMove and soluteMove
c###############################################################

      integer updatenrg,movnum,movtyp,isol,jsol,ranInt,savenrg
      double precision soluteEnergy,soluteEnergyD,extraEnergy,ranRange
      double precision ieold,ienew,thetaNew
      double precision inrgf,inrgb
      double precision thetaold,ranval,t1,t2,calcJPMF,corr
      double precision opmf,npmf,ocorr,ncorr
      character*4 statstring
      logical metropolis
 
      solThetaMove = .false.

c     Select a test particle at random
      isol = ranInt(1,NThetaSolutes)
      thetaold = SoluteTheta(isol)
 
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      
      NSolTheTry(isol) = NSolTheTry(isol) + 1

c     see whether we are updating the solute energies
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else 
        savenrg = DONTCARE
      endif

c     Check if a coupled solute must be moved as well
      jsol = DualTopologySolute(isol)
      if (jsol .lt. 0) jsol = -jsol

c     get the old energy of the solute - the internal angles were calculated
c     when it was built
      ieold = soluteEnergy(savenrg,isol,inrgb,inrgf)
     .            + extraEnergy(savenrg,inrgb,inrgf)

c     JM Dual Topology get energy of old solute if it is being moved
      if (jsol .ne. 0) ieold = ieold + soluteEnergy(savenrg,jsol,inrgb,inrgf)

c     Change the theta value of the molecule according to the random number
10    ranval = ranRange(-SolTDel(isol),SolTDel(isol))
c     Note that this uses ranRange, so it adds or removes +/- 0.15
      SoluteTheta(isol) = thetaold + ranval
      if (SoluteTheta(isol).gt.ONE) then
      goto 10
      endif
      if (SoluteTheta(isol).lt.ZERO) then 
      goto 10
      endif
      thetaNew = SoluteTheta(isol)

c     make the calculation go from old to new if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW

c     get the new energy
      ienew = soluteEnergy(savenrg,isol,inrgb,inrgf) + 
     .           extraEnergy(savenrg,inrgb,inrgf)

c     JM Dual Topology get energy of old solute if it is being moved
      if (jsol .ne. 0) ienew = ienew + soluteEnergy(savenrg,jsol,inrgb,inrgf)
    
c     Adds on corrections and solvation penalty
      t1 = 1.0-thetaold
      t2 = 1.0-SoluteTheta(isol)
      if (UseConcCorr) then
        corr = SolJCorr(isol)
      else
        corr = ZERO
      endif
      opmf = calcJPMF(isol,t1)*SoluteHydScale(isol)
      npmf = calcJPMF(isol,t2)*SoluteHydScale(isol)      
      ocorr = corr*thetaold
      ncorr = corr*SoluteTheta(isol)
      ieold = ieold -opmf + ocorr 
      ienew = ienew -npmf + ncorr

c     See if the new configurational energy needs to be updated...
      if (metropolis(ieold,ienew,movtyp)) then
        solThetaMove = .true.
        statstring = " 1 "          
        NSolTheAcp(isol) = NSolTheAcp(isol) + 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
           call updateSoluteEnergies(isol)
           call updateExtraEnergies
c          JM Dual Topology, update energies of second solute if necessary
           if (jsol .ne. 0) call updateSoluteEnergies(jsol)
        endif
      else
        solThetaMove = .false.
        statstring = " 0 "
        SoluteTheta(isol) = thetaold
      endif

!      if (printMoveEvery.eq.0) return
!      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,isol,SoluteName(isol)(1:10),ienew,
     &                        npmf,ncorr,ieold,opmf,ocorr,statstring,thetaNew,thetaold,ranval
12      format(i10,i2," ",a," ",6(g15.3,x),1x,a4,3f10.3)!" FTHETAMOVE "
        call printLine(MOVE,printstring)
cc       if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
!      endif
      
      return
      end
