      subroutine readChainTemplateBondLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c##################################################################
c
c     This routine reads in a chain bond template line
c
c     (C) Christopher Woods
c
c##################################################################

      character*(*) line
      integer itmpl
      character*4 at1,at2
      logical flex,dummy
      integer ipar0,ipar1,ibnd,findChainBond
      double precision delta
      
c     format is "bond at1 at2 (flex idelta) (dummy) (param ipar0 ipar1)"

      if (nwords.lt.3) then
        call printLine(WARNING,"Incorrectly formatted chain bond template line?")
        call printLine(WARNING,line)
      endif
      
c     read in the data
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)
      call sortTemplateBond(at1,at2)
      
      ibnd = findChainBond(at1,at2,itmpl)
      if (ibnd.ne.0) then
        call printLine(WARNING,"Overwriting chain template bond.")
        call printLine(WARNING,line)
      else
        ibnd = NChainBonds(itmpl) + 1
        if (ibnd.gt.MAXCHAINBONDS) then
          write(printstring,*) "Exceeded maximum number of chain template bonds (",
     .                          MAXCHAINBONDS,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXCHAINBONDS and recompile")
          return
        endif
        NChainBonds(itmpl) = ibnd
      endif
      
      ChainBndName(itmpl,ibnd,1) = at1
      ChainBndName(itmpl,ibnd,2) = at2
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      ChainBondType(itmpl,ibnd,1) = dummy
      ChainBondType(itmpl,ibnd,2) = flex
      ChainBondDelta(itmpl,ibnd) = delta
      ChainBondPar(itmpl,ibnd,1) = ipar0
      ChainBondPar(itmpl,ibnd,2) = ipar1
      
      return
      end
