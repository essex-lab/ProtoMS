      subroutine synchFEPSolutes(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'enums.inc'
c####################################################
c
c     This routine synchronises the coordinates
c     of solute isol with its forwards and
c     backwards copies
c
c     (C) Christopher Woods    14 August 2002
c
c####################################################

      integer isol,ifor,ibac,i,j
      double precision perturbForceField
      
c     make sure that we don't synch ignored solutes!
      if (IgnoreSolute(isol)) then
        write(printstring,*) "Why are we trying to synch with the ignored solute ",isol,"?"
        call printLine(WARNING,printstring)
        return
      endif
      
      ifor = NSolForwards(isol)
      ibac = NSolBackwards(isol)
      if (ifor.le.0 .or. ibac.le.0) then
        write(printstring,*) "Problem synching coordinates for solute ",isol,"(",
     .                       ifor,"-",ibac,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug or memory problem?")
        return
      endif
      
      do j=1,3
        do i=1,NSolAtoms(isol)
          ZSolVal(ifor,i,j) = ZSolVal(isol,i,j)
          ZSolVal(ibac,i,j) = ZSolVal(isol,i,j)
        enddo
        do i=1,3
          SolDummy(ifor,i,j) = SolDummy(isol,i,j)
          SolDummy(ibac,i,j) = SolDummy(isol,i,j)
        enddo
      enddo

c     reapply the variable geometries
      do i=1,NSolVariable(isol)
        ZSolVal(ifor,SolVariableAtm(isol,i),SolVariableTyp(isol,i)) = 
     .   perturbForceField(SolVariableVal(isol,i,1),SolVariableVal(isol,i,2),lambdaF)
        ZSolVal(ibac,SolVariableAtm(isol,i),SolVariableTyp(isol,i)) = 
     .   perturbForceField(SolVariableVal(isol,i,1),SolVariableVal(isol,i,2),lambdaB)
        ZSolVal(isol,SolVariableAtm(isol,i),SolVariableTyp(isol,i)) = 
     .   perturbForceField(SolVariableVal(isol,i,1),SolVariableVal(isol,i,2),lambda)        
      enddo
     
c     rebuild the solutes - this recalculates all of the internal angles
      call buildSolute(isol)
      call buildSolute(ifor)
      call buildSolute(ibac)
      
      return
      end
