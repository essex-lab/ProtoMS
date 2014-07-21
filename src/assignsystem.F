      subroutine assignSystem
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'pressure.inc'
c#################################################################
c
c       This routine controls the assignment of templates
c       for the loaded system and builds all of the data
c       structures
c
c       (C) Christopher Woods, October 2004
c
c#################################################################

      integer i
      logical ok
      
c     assign the solutes first as they are quickest and most likely 
c     to contain user errors (thus the user needs quick feedback of a
c     problem)
      if (NSolutes.ne.0) then
c       add solute dummies and assign molecule templates
        do i=1,NSolutes
c         assign the molecule templates
          call assignSoluteTemplates(i)
c         now get the centre and radius of this solute
          call getSoluteCentreRadius(i)
        enddo
      endif
        
      if (NGCsolutes.gt.0) then
c       now assign the templates of all of the solvents
        call assignGCsoluteTemplates

c       now all of the solvents have been loaded, get their centres of
c       geometry, and internal axes
c        write(printstring,*) "Calculating solvent axes..."
c        call printLine(INFO,printstring)
        do i=1,NGCsolutes
          call addGCsoluteAxes(i)
        enddo
      endif


       
c     assign and build the zmatricies for all of the proteins      
      if (NProteins.gt.0) then      
        do i=1,NProteins
          call buildProtein(i)
        enddo
      endif

c     sort out the boundary conditions for the solvent and assign templates

c     if there are no solvents then there are no boundaries!
      if (NSolvents.le.0) then
        SolventBox = .false.
        SolventCap = .false.
        boundary = NOBOUNDARIES
        call prettyLine(INFO,"Using no boundary conditions (gas-phase)")
      else if (SolventBox) then
        if (boundary.eq.AGREESOLVENT) boundary = PERIODICBOUNDARIES
        
        if (boundary.eq.PERIODICBOUNDARIES) then
          call prettyLine(INFO,"Using periodic boundary conditions")
        else if (boundary.eq.CUBICHARMONIC) then
          call prettyLine(INFO,"Using cubic harmonic boundary conditions")              
        endif            

        ok = .true.
        do i=1,3
          boxDimension(i) = boxTop(i)-boxOrigin(i)
          boxDimension2(i) = HALF*boxDimension(i)
          if (boxDimension(i).le.ZERO) ok = .false.
        enddo
        
        if (.not.ok) then
          call printLine(WARNING,"Serious problem with box dimensions!")
          write(printstring,20) (BoxOrigin(i),i=1,3)
          call printLine(INFO,printstring)
          write(printstring,21) (BoxDimension(i),i=1,3)
          call printLine(INFO,printstring)
          write(printstring,22) (BoxTop(i),i=1,3)
          call printLine(INFO,printstring)
          call printLine(WARNING,"Will recalculate boundaries - these will need to be equilibrated!")
          call calculateSolventBox
        endif
        
        write(printstring,20) (BoxOrigin(i),i=1,3)
        call printLine(INFO,printstring)
        write(printstring,21) (BoxDimension(i),i=1,3)
        call printLine(INFO,printstring)
        write(printstring,22) (BoxTop(i),i=1,3)
        call printLine(INFO,printstring)
20      format("Solvent box origin = (",f8.3,",",f8.3,",",f8.3,")")
21      format("Solvent box dimensions = (",f8.3,",",f8.3,",",f8.3,")")
22      format("Solvent box top = (",f8.3,",",f8.3,",",f8.3,")")
      
c       make sure that the box is greater than the cutoff
        ok = .true.
        do i=1,3
          if (cut.ge.HALF*BoxDimension(i)) ok = .false.
        enddo
        
        if (.not.ok) then
          write(printstring,*) "Non-bonded cutoff (",cut," A) is too large for this solvent box!"
          call prettyLine(INFO,printstring)
          call printLine(WARNING,printstring)
        endif
        
c       calculate the volume of the box
        volume = BoxDimension(1)*BoxDimension(2)*BoxDimension(3)
        write(printstring,*) "Initial volume = ",volume," A^3"
        call prettyLine(INFO,printstring)
      else if (SolventCap) then
        if (boundary.eq.AGREESOLVENT) boundary = SPHERICALHARMONIC
        
        call prettyLine(INFO,"Using solvent cap boundary conditions")
        write(printstring,30) (cap(i),i=1,3)
        call printLine(INFO,printstring)
        write(printstring,31) capradius,capconstant
        call printLine(INFO,printstring)
30      format('Solvent cap origin = (',f8.3,',',f8.3,',',f8.3,')')
31      format('Solvent cap radius = ',f8.3,' A, constant = ',f8.3,' kcal mol-1 A-2')

        if (capradius.le.ZERO .or. capconstant.le.ZERO) then
          call prettyLine(WARNING,"Are you sure that the parameters of the solvent cap are correct?")
        endif

      else if (boundary.eq.PERIODICBOUNDARIES .or. boundary.eq.CUBICHARMONIC) then
        if (boundary.eq.PERIODICBOUNDARIES) then
          call prettyLine(INFO,"Using periodic boundary conditions")
        else if (boundary.eq.CUBICHARMONIC) then
          call prettyLine(INFO,"Using cubic harmonic boundary conditions")              
        endif            
        call printLine(INFO,"Automatically generating solvent box dimensions - you will need to")
        call printLine(INFO,"run some NPT equilibration to get these fully right")                

        call calculateSolventBox

        write(printstring,20) (BoxOrigin(i),i=1,3)
        call printLine(INFO,printstring)
        write(printstring,21) (BoxDimension(i),i=1,3)
        call printLine(INFO,printstring)
        write(printstring,22) (BoxTop(i),i=1,3)
        call printLine(INFO,printstring)
      
c       calculate the volume of the box
        volume = BoxDimension(1)*BoxDimension(2)*BoxDimension(3)
        write(printstring,*) "Initial volume = ",volume," A^3"
        call prettyLine(INFO,printstring)
      else if (boundary.eq.SPHERICALHARMONIC) then
        call printLine(WARNING,"Requested solvent cap boundaries but have not provided any information")
        call printLine(WARNING,"to describe the solvent cap! Turning solvent cap boundaries off!")
        boundary = NOBOUNDARIES
      else if (boundary.eq.AGREESOLVENT) then
        call printLine(WARNING,"Requested boundaries that agree with the solvent, but the solvent doesn't")
        call printLine(WARNING,"specify any sort of boundaries! Using gas-phase boundaries...")
        boundary = NOBOUNDARIES
      endif

      if (NSolvents.gt.0) then
       
c       now assign the templates of all of the solvents
        call assignSolventTemplates

c       now all of the solvents have been loaded, get their centres of
c       geometry, and internal axes
        write(printstring,*) "Calculating solvent axes..."
        call printLine(INFO,printstring)
        do i=1,NSolvents
          call addSolventAxes(i)
        enddo

        write(printstring,*) "All axes successfully generated!"
        call printLine(INFO,printstring)
      endif
        
      end
