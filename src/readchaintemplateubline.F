      subroutine readChainTemplateUBLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
c###########################################################
c
c     This routine reads in a chain Urey-Bradley template line
c
c     (C) Christopher Woods
c
c###########################################################

      character*(*) line
      integer itmpl
      character*4 at1,at2,at3
      logical flex,dummy
      integer ipar0,ipar1,iuby,findChainUreyBradley
      double precision delta
      
c     format is "ureybradley atm1 atm2 atm3 (dummy) (param ipar0 ipar1)
c      (ubys cannot be flexible)

      if (nwords.lt.4) then
        call printLine(WARNING,"Incorrectly formatted chain Urey-Bradley template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)
      at3 = words(4)(1:4)

      call sortTemplateUreyBradley(at1,at2,at3)
      iuby = findChainUreyBradley(at1,at2,at3,itmpl)
      if (iuby.eq.0) then
        iuby = NChainUreyBradleys(itmpl) + 1
        if (iuby.gt.MAXCHAINUREYBRADLEYS) then
          write(printstring,*) "Exceeded maximum number of chain template Urey-Bradley terms (",
     .                         MAXCHAINUREYBRADLEYS,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXCHAINUREYBRADLEYS and recompile")
          return
        endif
        NChainUreyBradleys(itmpl) = iuby
      else
        write(printstring,*) "Overwriting chain template Urey-Bradley for ",at1,"-",
     .                       at2,"-",at3
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
            
      ChainUByName(itmpl,iuby,1) = at1
      ChainUByName(itmpl,iuby,2) = at2
      ChainUByName(itmpl,iuby,3) = at3
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      ChainUreyBradleyType(itmpl,iuby) = dummy
      if (flex) then
        call printLine(WARNING,"Cannot have flexible Urey-Bradleys!")
        call printLine(WARNING,line)
      endif
      
      ChainUreyBradleyPar(itmpl,iuby,1) = ipar0
      ChainUreyBradleyPar(itmpl,iuby,2) = ipar1
      
      return
      end  
