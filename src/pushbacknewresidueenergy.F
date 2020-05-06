      subroutine pushBackNewResidueEnergy
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'stacks.inc'
c###########################################################
c
c     This subroutine pushes back the new residue energies
c     onto the stack
c
c     (C) Christopher Woods   8 February 2002
c
c###########################################################

      integer i
      
      newStackDepth = newStackDepth + 1
      if (newStackDepth.gt.MAXRESSTACKDEPTH) then
        write(printstring,*) "Exceeded residue energy stack depth! (",MAXRESSTACKDEPTH,
     .                       ") - increase MAXRESSTACKDEPTH and recompile"
        call closeProgram(-1,printstring)
      endif
      
c     save the newresidue energies
      if (NProteins.gt.1) then
        do i=1,NProteins
          newStackResProC(newStackDepth,i) = newResProC(i)
          newStackResProL(newStackDepth,i) = newResProL(i)
          newStackResProCF(newStackDepth,i) = newResProCF(i)
          newStackResProLF(newStackDepth,i) = newResProLF(i)
          newStackResProCB(newStackDepth,i) = newResProCB(i)
          newStackResProLB(newStackDepth,i) = newResProLB(i)
          newStackResProDist(newStackDepth,i) = newResProDist(i)
        enddo
      endif
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          newStackResSolC(newStackDepth,i) = newResSolC(i)
          newStackResSolL(newStackDepth,i) = newResSolL(i)
          newStackResSolCF(newStackDepth,i) = newResSolCF(i)
          newStackResSolLF(newStackDepth,i) = newResSolLF(i)
          newStackResSolCB(newStackDepth,i) = newResSolCB(i)
          newStackResSolLB(newStackDepth,i) = newResSolLB(i)
          newStackResSolDist(newStackDepth,i) = newResSolDist(i)
        enddo
      endif
      newStackResSvnC(newStackDepth) = newResSvnC
      newStackResSvnL(newStackDepth) = newResSvnL
      newStackResSvnCF(newStackDepth) = newResSvnCF
      newStackResSvnLF(newStackDepth) = newResSvnLF
      newStackResSvnCB(newStackDepth) = newResSvnCB
      newStackResSvnLB(newStackDepth) = newResSvnLB
      newStackResSvnDist(newStackDepth) = newResSvnDist
      
      newStackResBnd(newStackDepth) = newResBnd
      newStackResAng(newStackDepth) = newResAng
      newStackResDih(newStackDepth) = newResDih
      newStackResUBy(newStackDepth) = newResUBy
      newStackResNBC(newStackDepth) = newResNBC
      newStackResNBL(newStackDepth) = newResNBLJ
      
      newStackResBndF(newStackDepth) = newResBndF
      newStackResBndB(newStackDepth) = newResBndB
      newStackResAngF(newStackDepth) = newResAngF
      newStackResAngB(newStackDepth) = newResAngB
      newStackResDihF(newStackDepth) = newResDihF
      newStackResDihB(newStackDepth) = newResDihB
      newStackResUByF(newStackDepth) = newResUByF
      newStackResUByB(newStackDepth) = newResUByB
      newStackResNBCF(newStackDepth) = newResNBCF
      newStackResNBCB(newStackDepth) = newResNBCB
      newStackResNBLF(newStackDepth) = newResNBLJF
      newStackResNBLB(newStackDepth) = newResNBLJB
      
      newStackResNrg(newStackDepth) = newResNrg  
      newStackResNrgF(newStackDepth) = newResNrgF
      newStackResNrgB(newStackDepth) = newResNrgB
      
      return
      end
