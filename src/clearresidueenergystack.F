      subroutine clearResidueEnergyStack
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'stacks.inc'
c###################################################################
c
c     This routine clears the residue energy stack
c      (shallow clear - it only sets the counters back to 0)
c
c     (C) Christopher Woods   8 February 2002
c
c###################################################################

      if (newStackDepth.ne.oldStackDepth) then
        write(printstring,*) "The residue energy stacks are not aligned! ",newStackDepth,oldStackDepth
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is a severe programming bug!")
        call closeProgram(-1,"ProtoMS is exiting because of a severe programming bug!")
      endif

      newStackDepth = 0
      oldStackDepth = 0
      
      return
      end
