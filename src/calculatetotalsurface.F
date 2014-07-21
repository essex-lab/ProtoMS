      subroutine calculateTotalSurface(showSurface)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'surface.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'enums.inc'
c################################################
c     Calculate the Surface Area of the system 
c     from scratch
c     and prints a detailed output if requested
c     
c     Julien Michel November 2004
c################################################
      logical showSurface
      integer i,j,nres
      double precision TotalSurface
      double precision surface,surfacef,surfaceb
      double precision sanrg,sanrgf,sanrgb

      call clearSurfaces()
      
      call setAllMoved()

      surface = Totalsurface(NEW,surfacef,surfaceb,sanrg,sanrgf,sanrgb)

      call clearMoveList()

      if (showSurface) then 
        if (NProteins .gt. ZERO) then      
          do i=1,NProteins
            nres = NResidues(i)
            do j=1,nres
              write(printstring,31) 'Protein ',i,' Residue ',j,NewResidueSurf(i,j),NewResidueSurfF(i,j),NewResidueSurfB(i,j)
              call printLine(INFO,printstring)
            enddo
          enddo
        endif
        
        if (NSolutes .gt. ZERO) then
          do i=1,NSolutes
            if (IgnoreSolute(i)) continue
            write(printstring,32) 'Solute ',i,NewSoluteSurf(i),NewSoluteSurfF(i),NewSoluteSurfB(i)
            call printLine(INFO,printstring)
          enddo
        endif
        
        if (NSolvents .gt. ZERO) then
          do i=1,NSolvents
            write(printstring,32) 'Solvent ',i,NewSolventSurf(i),NewSolventSurfF(i),NewSolventSurfB(i)
            call printLine(INFO,printstring)
          enddo
        endif

 31     format(a10,i3,a10,i3,' Surface (Forwards,Backwards) : ',f8.2,'(',f8.2,')','(',f8.2,')')
 32     format(a10,i3,' Surface (Forwards,Backwards) : ',f8.2,'(',f8.2,')','(',f8.2,')')
        
        write(printstring,33) surface,surfacef,surfaceb
 33     format('Surface        (Forwards,Backwards) : ',f8.2,'(',f8.2,')','(',f8.2,')')
        call printLine(INFO,printstring)
        write(printstring,34) sanrg,sanrgf,sanrgb
 34     format('Tension Energy (Forwards,Backwards) : ',f8.3,'(',f8.3,')','(',f8.3,')')
        call printLine(INFO,printstring)
      endif

      return
      end
