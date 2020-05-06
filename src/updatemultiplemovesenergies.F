      subroutine updateMultipleMovesEnergies
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'movelist.inc'
      include 'stacks.inc'
      include 'printer.inc'
c########################################
c
c     This routine updates the energy 
c     components of all the parts 
c     of the system that have been moved 
c     (typically after a multiple move)
c
c     (c) Julien Michel May 2005
c########################################
      
      integer i,moved
      logical hasSoluteMoved,hasSolventMoved

c$$$      if (NProteins .gt. 0) then
c$$$        do 10 i=1,NProteins
c$$$          nres = NResidues(i)
c$$$          do 12 j=1,nres
c$$$            if ( hasResidueMoved(i,j) ) then 
c$$$c             Humm, how do I get the energy components now !!
c$$$c             It looks like I need to save them somewhere !!!
c$$$c             The stack is not appropriate here as I will not 
c$$$c             necessarily acces the residues in the right order
c$$$c             (even though it may actually work because the length
c$$$c             of the stack should be equal to the number of moved 
c$$$c             residues)
c$$$c             so I will instead, repeatidly pop the stack until it is 
c$$$c             clear. It is also a bit faster
c$$$c             It means that I should clear the residueenergy stack
c$$$c             when I am clearing the movelist
c$$$              call popAllResidueEnergy
c$$$              call updateResidueEnergies(i)
c$$$            endif              
c$$$ 12       continue
c$$$ 10     continue
c$$$      endif
      moved = newStackDepth
c      write(printstring,*) "NUMBER OF MOVED RESIDUES ",moved
c      call printLine(DEBUG,printstring)
      do  i=1,moved
        call popAllResidueEnergy
c       Problem, when I pop a residue, I must know from which protein 
c       that residue came from 
        call updateResidueEnergies(1)
      enddo
      
      if (NSolutes.gt. 0) then
        do 20 i=1,NSolutes
          if ( hasSoluteMoved(i) ) then
c           Where do I get the solute energies components ?
            call updateSoluteEnergies(i)
          endif
 20     continue
      endif

      if (NSolvents .gt. 0) then
        do 40 i=1,NSolvents
          if ( hasSolventMoved(i) ) then
c           I don't think I need to have saved all the components...
            call updateSolventEnergies()
c           ...and I should do it once 
            goto 41
          endif
 40     continue
 41     continue
      endif      

      call updateextraenergies()


      return
      end
