      subroutine setInitialTitrationState
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#################################################################
c
c       This routine is called to set the initial titration
c       state of the proteins (if 'neutralise' is set)
c
c       (C) Christopher Woods   9 August 2002
c
c#################################################################

      logical protonate,deprotonate
      integer i,ichg,ires
      
c     if there are any proteins...
      if (NProteins.gt.0) then
c       if requested, neutralise the protein...
        if (NeutraliseProtein) then
          write(printstring,*) "Neutralising the proteins..."
          call printLine(INFO,printstring)
c         loop over all of the proteins, and neutralise them...        
          do 100 i=1,NProteins
            ichg = Int(TotalCharge(i))
            write(printstring,*) "Total charge on protein ",i," is ",TotalCharge(i)
            call printLine(INFO,printstring)
            if (ichg.lt.desiredCharge) then
c             we need to add protons...
101           continue
              if (NBasicRes(i).gt.0) then
c               protonate the first few basic residues...
                ires = BasicRes(i,1)
                if (protonate(i,ires)) then
                   write(printstring,10) ResNam(i,ires),ResNum(i,ires),i
10                 format("Protonating ",a,"(",i3,") of protein ",i2)
                   call printLine(DETAIL,printstring)
                   TotalCharge(i) = TotalCharge(i) + 1
                   ichg = ichg + 1
                else
                  write(printstring,*) "Problem protonating ",ResNam(i,ires),ResNum(i,ires),i
                  call closeProgram(-1,printstring)
                endif
                if (TotalCharge(i).lt.desiredCharge) goto 101
              endif                
            else if (ichg.gt.desiredCharge) then
c             we need to take protons...
102           continue
              if (NAcidicRes(i).gt.0) then
c               deprotonate the first few acidic residues...
                ires = AcidicRes(i,1)
                if (deprotonate(i,ires)) then
                  write(printstring,11) ResNam(i,ires),ResNum(i,ires),i
11                 format("Deprotonating ",a,"(",i3,") of protein ",i2)
                   call printLine(DETAIL,printstring)
                   TotalCharge(i) = TotalCharge(i) - 1
                   ichg = ichg - 1
                else
                  write(printstring,*) "Problem deprotonating ",ResNam(i,ires),ResNum(i,ires),i
                  call closeProgram(-1,printstring)
                endif
                if (TotalCharge(i).gt.desiredCharge) goto 102
              endif
            endif  
c           resum the charge of this protein...
            call sumCharge(i)
c           write out some info...
            write (printstring,12) NAcidicRes(i),NBasicRes(i)
12          format("There are now ",i4," acidic residues, and ",i4," basic residues...")
            call printLine(INFO,printstring)
100       continue            
        endif
      endif
                                   
      end
