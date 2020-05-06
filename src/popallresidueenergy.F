      subroutine popAllResidueEnergy
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'stacks.inc'
c##########################################################
c
c     This subroutine pops the old and new residue
c     energies off of the stack (LIFO)
c
c     (C) Christopher Woods   8 February 2002
c
c##########################################################

      integer i
      
      if (newStackDepth.ne.oldStackDepth) then
        write(printstring,*) "The old and new residue energy stacks are",
     .      " not aligned! ",newStackDepth,oldStackDepth
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is a severe programming bug!")
        call closeProgram(-1,"ProtoMS is exiting because of a severe programming bug!")
      else if (newStackDepth.le.0) then
        write(printstring,*) "We have popped off the end of the residue energy stacks!"
        call printLine(WARNING,"This is a severe programming bug!")
        call closeProgram(-1,"ProtoMS is exiting because of a severe programming bug!")
      endif
      
c     first extract the new and old energies...
      if (NProteins.gt.1) then
        do i=1,NProteins
          newResProC(i) = newStackResProC(newStackDepth,i)
          newResProL(i) = newStackResProL(newStackDepth,i)
          newResProCF(i) = newStackResProCF(newStackDepth,i)
          newResProLF(i) = newStackResProLF(newStackDepth,i)
          newResProCB(i) = newStackResProCB(newStackDepth,i)
          newResProLB(i) = newStackResProLB(newStackDepth,i)
          newResProDist(i) = newStackResProDist(newStackDepth,i)
          oldResProC(i) = oldStackResProC(newStackDepth,i)
          oldResProL(i) = oldStackResProL(newStackDepth,i)
          oldResProCF(i) = oldStackResProCF(newStackDepth,i)
          oldResProLF(i) = oldStackResProLF(newStackDepth,i)
          oldResProCB(i) = oldStackResProCB(newStackDepth,i)
          oldResProLB(i) = oldStackResProLB(newStackDepth,i)
          oldResProDist(i) = oldStackResProDist(newStackDepth,i)
        enddo
      endif
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          newResSolC(i) = newStackResSolC(newStackDepth,i)
          newResSolL(i) = newStackResSolL(newStackDepth,i)
          newResSolCF(i) = newStackResSolCF(newStackDepth,i)
          newResSolLF(i) = newStackResSolLF(newStackDepth,i)
          newResSolCB(i) = newStackResSolCB(newStackDepth,i)
          newResSolLB(i) = newStackResSolLB(newStackDepth,i)
          newResSolDist(i) = newStackResSolDist(newStackDepth,i)
          oldResSolC(i) = oldStackResSolC(newStackDepth,i)
          oldResSolL(i) = oldStackResSolL(newStackDepth,i)
          oldResSolCF(i) = oldStackResSolCF(newStackDepth,i)
          oldResSolLF(i) = oldStackResSolLF(newStackDepth,i)
          oldResSolCB(i) = oldStackResSolCB(newStackDepth,i)
          oldResSolLB(i) = oldStackResSolLB(newStackDepth,i)
          oldResSolDist(i) = oldStackResSolDist(newStackDepth,i)
        enddo
      endif
      
      newResSvnC = newStackResSvnC(newStackDepth)
      newResSvnL = newStackResSvnL(newStackDepth)
      newResSvnCF = newStackResSvnCF(newStackDepth)
      newResSvnLF = newStackResSvnLF(newStackDepth)
      newResSvnCB = newStackResSvnCB(newStackDepth)
      newResSvnLB = newStackResSvnLB(newStackDepth)
      newResSvnDist = newStackResSvnDist(newStackDepth)
      oldResSvnC = oldStackResSvnC(newStackDepth)
      oldResSvnL = oldStackResSvnL(newStackDepth)
      oldResSvnCF = oldStackResSvnCF(newStackDepth)
      oldResSvnLF = oldStackResSvnLF(newStackDepth)
      oldResSvnCB = oldStackResSvnCB(newStackDepth)
      oldResSvnLB = oldStackResSvnLB(newStackDepth)
      oldResSvnDist = oldStackResSvnDist(newStackDepth)
      
      newResBnd = newStackResBnd(newStackDepth)
      newResAng = newStackResAng(newStackDepth)
      newResDih = newStackResDih(newStackDepth)
      newResUBy = newStackResUBy(newStackDepth)
      newResNBC = newStackResNBC(newStackDepth)
      newResNBLJ = newStackResNBL(newStackDepth)
      
      oldResBnd = oldStackResBnd(newStackDepth)
      oldResAng = oldStackResAng(newStackDepth)
      oldResDih = oldStackResDih(newStackDepth)
      oldResUBy = oldStackResUBy(newStackDepth)
      oldResNBC = oldStackResNBC(newStackDepth)
      oldResNBLJ = oldStackResNBL(newStackDepth)
      
      newResNrg = newStackResNrg(newStackDepth)
      oldResNrg = oldStackResNrg(newStackDepth)  
      
      newResNrgF = newStackResNrgF(newStackDepth)
      oldResNrgF = oldStackResNrgF(newStackDepth)
      newResNrgB = newStackResNrgB(newStackDepth)
      oldResNrgB = oldStackResNrgB(newStackDepth)

      newResBndF = newStackResBndF(newStackDepth)
      newResBndB = newStackResBndB(newStackDepth)
      newResAngF = newStackResAngF(newStackDepth)
      newResAngB = newStackResAngB(newStackDepth)
      newResDihF = newStackResDihF(newStackDepth)
      newResDihB = newStackResDihB(newStackDepth)
      newResUByF = newStackResUByF(newStackDepth)
      newResUByB = newStackResUByB(newStackDepth)
      newResNBCF = newStackResNBCF(newStackDepth)
      newResNBCB = newStackResNBCB(newStackDepth)
      newResNBLJF = newStackResNBLF(newStackDepth)
      newResNBLJB = newStackResNBLB(newStackDepth)

      oldResBndF = oldStackResBndF(oldStackDepth)
      oldResBndB = oldStackResBndB(oldStackDepth)
      oldResAngF = oldStackResAngF(oldStackDepth)
      oldResAngB = oldStackResAngB(oldStackDepth)
      oldResDihF = oldStackResDihF(oldStackDepth)
      oldResDihB = oldStackResDihB(oldStackDepth)
      oldResUByF = oldStackResUByF(oldStackDepth)
      oldResUByB = oldStackResUByB(oldStackDepth)
      oldResNBCF = oldStackResNBCF(oldStackDepth)
      oldResNBCB = oldStackResNBCB(oldStackDepth)
      oldResNBLJF = oldStackResNBLF(oldStackDepth)
      oldResNBLJB = oldStackResNBLB(oldStackDepth)

      newStackDepth = newStackDepth - 1
      oldStackDepth = oldStackDepth - 1
      
      return
      end
