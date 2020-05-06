      double precision function lambdaderivative()
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'movelist.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'quickpot.inc'
c###########################################################
c
c     This function returns the derivative of the total energy
c     with respect to lambda
c
c     Only works with dual-topology
c
c     (C) Samuel Genheden 30 July 2014
c
c###########################################################

      integer savenrg
      integer i,j
      double precision nrg,nrgb,nrgf
      double precision cnrg,ljnrg,mindist
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision nbcnrg,nbljnrg
      double precision nbcnrgf,nbljnrgf
      double precision nbcnrgb,nbljnrgb
      double precision scl,sclf,sclb
      double precision getLowScale,getHighScale
      double precision soluteSoluteEnergy
      double precision proteinSoluteEnergy
      double precision soluteSolventEnergy

      lambdaderivative = ZERO
      
      doDerivative = .true.
      CoulSoftPower = CoulSoftPower + 1.0
c     now do the same for all of the solutes...

      if (NSolutes.gt.1) then
        do i=1,NSolutes-1
          if (.not.IgnoreSolute(i)) then
            do j=i+1,NSolutes
              if (.not.IgnoreSolute(j)) then
                nrg = soluteSoluteEnergy(i,j,
     .                      cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .                      scl,sclf,sclb)
                lambdaderivative = lambdaderivative + nrg              
              endif
            enddo
          endif
        enddo
      endif
      if (NSolutes.gt.0) then
        do j=1,NSolutes
          if ((.not.IgnoreSolute(j)).and.(DualTopologySolute(j).ne.0)) then
c           see if this is a dual topology solute
            if (DualTopologySolute(j).gt.0) then
              scl = getLowScale(Lambda)
              sclb = getLowScale(LambdaB)
              sclf = getLowScale(LambdaF)
              doLowDerivative = .True.
            else if (DualTopologySolute(j).eq.0) then
              scl = ONE
              sclb = ONE
              sclf = ONE
            else
              scl = getHighScale(Lambda)
              sclb = getHighScale(LambdaB)
              sclf = getHighScale(LambdaF)
              doLowDerivative = .False.
            endif
          
c           solute-protein energies...
            if (NProteins.gt.0) then
              do i=1,NProteins
                nrg = proteinSoluteEnergy(DONTCARE,i,j,
     .                      cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .                      scl,sclf,sclb)
                lambdaderivative = lambdaderivative + nrg
              enddo
            endif
c           solute-solvent energies...
            nrg = soluteSolventEnergy(DONTCARE,j,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .              scl,sclf,sclb)
            lambdaderivative = lambdaderivative + nrg
          endif
        enddo
      endif

      doDerivative = .false.
      CoulSoftPower = CoulSoftPower - 1.0
      return

      end
