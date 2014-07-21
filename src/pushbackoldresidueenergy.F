      subroutine pushBackOldResidueEnergy
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'stacks.inc'
c###########################################################
c
c     This subroutine pushes back the old residue energies
c     onto the stack
c
c     (C) Christopher Woods   8 February 2002
c
c###########################################################

      integer i

      oldStackDepth = oldStackDepth + 1
      if (oldStackDepth.gt.MAXRESSTACKDEPTH) then
        write(printstring,*) "Exceeded residue energy stack depth! (",MAXRESSTACKDEPTH,
     .                       ") - increase MAXRESSTACKDEPTH and recompile"
        call closeProgram(-1,printstring)
      endif
      
c     now save the oldresidue energies
      if (NProteins.gt.1) then
        do i=1,NProteins
          oldStackResProC(oldStackDepth,i) = oldResProC(i)
          oldStackResProL(oldStackDepth,i) = oldResProL(i)
          oldStackResProCF(oldStackDepth,i) = oldResProCF(i)
          oldStackResProLF(oldStackDepth,i) = oldResProLF(i)
          oldStackResProCB(oldStackDepth,i) = oldResProCB(i)
          oldStackResProLB(oldStackDepth,i) = oldResProLB(i)
          oldStackResProDist(oldStackDepth,i) = oldResProDist(i)
        enddo
      endif
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          oldStackResSolC(oldStackDepth,i) = oldResSolC(i)
          oldStackResSolL(oldStackDepth,i) = oldResSolL(i)
          oldStackResSolCF(oldStackDepth,i) = oldResSolCF(i)
          oldStackResSolLF(oldStackDepth,i) = oldResSolLF(i)
          oldStackResSolCB(oldStackDepth,i) = oldResSolCB(i)
          oldStackResSolLB(oldStackDepth,i) = oldResSolLB(i)
          oldStackResSolDist(oldStackDepth,i) = oldResSolDist(i)
        enddo
      endif
      oldStackResSvnC(oldStackDepth) = oldResSvnC
      oldStackResSvnL(oldStackDepth) = oldResSvnL
      oldStackResSvnCF(oldStackDepth) = oldResSvnCF
      oldStackResSvnLF(oldStackDepth) = oldResSvnLF
      oldStackResSvnCB(oldStackDepth) = oldResSvnCB
      oldStackResSvnLB(oldStackDepth) = oldResSvnLB
      oldStackResSvnDist(oldStackDepth) = oldResSvnDist
      
      oldStackResBnd(oldStackDepth) = oldResBnd
      oldStackResAng(oldStackDepth) = oldResAng
      oldStackResDih(oldStackDepth) = oldResDih
      oldStackResUBy(oldStackDepth) = oldResUBy
      oldStackResNBC(oldStackDepth) = oldResNBC
      oldStackResNBL(oldStackDepth) = oldResNBLJ
      
      oldStackResBndF(oldStackDepth) = oldResBndF
      oldStackResBndB(oldStackDepth) = oldResBndB
      oldStackResAngF(oldStackDepth) = oldResAngF
      oldStackResAngB(oldStackDepth) = oldResAngB
      oldStackResDihF(oldStackDepth) = oldResDihF
      oldStackResDihB(oldStackDepth) = oldResDihB
      oldStackResUByF(oldStackDepth) = oldResUByF
      oldStackResUByB(oldStackDepth) = oldResUByB
      oldStackResNBCF(oldStackDepth) = oldResNBCF
      oldStackResNBCB(oldStackDepth) = oldResNBCB
      oldStackResNBLF(oldStackDepth) = oldResNBLJF
      oldStackResNBLB(oldStackDepth) = oldResNBLJB

      oldStackResNrg(oldStackDepth) = oldResNrg
      oldStackResNrgF(oldStackDepth) = oldResNrgF
      oldStackResNrgB(oldStackDepth) = oldResNrgB
      
      return
      end
