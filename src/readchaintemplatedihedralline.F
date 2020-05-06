      subroutine readChainTemplateDihedralLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c########################################################
c
c     This routine reads in a chain dihedral template line
c
c     (C) Christopher Woods
c
c########################################################

      character*(*) line
      integer itmpl
      character*4 at1,at2,at3,at4
      logical flex,dummy
      integer ipar0,ipar1,idih,findChainDihedral
      double precision delta,deg2Rad
      
c     format is "dihedral atm1 atm2 atm3 atm4 (flex idelta) (dummy) (param ipar0 ipar1)

      if (nwords.lt.5) then
        call printLine(WARNING,"Incorrectly formatted chain dihedral template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)
      at3 = words(4)(1:4)
      at4 = words(5)(1:4)
            
      call sortTemplateDihedral(at1,at2,at3,at4)
      idih = findChainDihedral(at1,at2,at3,at4,itmpl)
      if (idih.eq.0) then
        idih = NChainDihedrals(itmpl) + 1
        if (idih.gt.MAXCHAINDIHEDRALS) then
          write(printstring,*) "Exceeded maximum number of chain template dihedrals (",
     .                         MAXCHAINDIHEDRALS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXCHAINDIHEDRALS and recompile"
          call closeProgram(-1,printstring)
          return
        endif 
        NChainDihedrals(itmpl) = idih
      else
        write(printstring,*) "Overwriting chain template angle for ",at1,"-",at2,
     .                      "-",at3,"-",at4
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
      
      ChainDihName(itmpl,idih,1) = at1
      ChainDihName(itmpl,idih,2) = at2
      ChainDihName(itmpl,idih,3) = at3
      ChainDihName(itmpl,idih,4) = at4
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      ChainDihedralType(itmpl,idih,1) = dummy
      ChainDihedralType(itmpl,idih,2) = flex
      ChainDihedralDelta(itmpl,idih) = deg2Rad(delta)
      ChainDihedralPar(itmpl,idih,1) = ipar0
      ChainDihedralPar(itmpl,idih,2) = ipar1
      
      return
      end
