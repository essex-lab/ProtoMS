      subroutine readSolTemplateDihedralLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###############################################################
c
c     This routine reads in the dihedrals from the molecule template
c
c     (C) Christopher Woods 
c
c###############################################################

      character*(*) line
      integer itmpl
      character*4 iatm1,ires1,iatm2,ires2,iatm3,ires3,iatm4,ires4
      logical dummy,flex
      double precision delta
      integer ipar0,ipar1
      integer findMolTemplateDihedral
      integer idih
      double precision deg2Rad
      double precision well1,well2,well3,well4
      integer nwells
c     format is;
c     'dihedral iatm1 ires1 iatm2 ires2 iatm3 ires3 iatm4 ires4 (dummy) (flex delta) (param ipar0 ipar1)'
      
      if (nwords.lt.9) then
        call printLine(WARNING,"Incorrectly formatted solute template dihedral line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the data
      iatm1 = words(2)(1:4)
      ires1 = words(3)(1:4)
      iatm2 = words(4)(1:4)
      ires2 = words(5)(1:4)
      iatm3 = words(6)(1:4)
      ires3 = words(7)(1:4)
      iatm4 = words(8)(1:4)
      ires4 = words(9)(1:4)
            
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
c     JM mod, check if some wells have been defined
      call getDihWells(well1,well2,well3,well4,nwells)
      

c     see if the dihedral already exists
      idih = findMolTemplateDihedral(iatm1,ires1,iatm2,ires2,
     .                               iatm3,ires3,iatm4,ires4,itmpl)
      if (idih.ne.0) then
        write(printstring,*) "Overwriting information for solute dihedral (",
     .                       iatm1,",",ires1,")-(",iatm2,",",ires2,")-(",
     .                       iatm3,",",ires4,")-(",iatm4,",",ires4,")"
        call printLine(INFO,printstring)
        call printLine(INFO,line)
      else
        idih = NMolTemplateDihedrals(itmpl) + 1
        if (idih.gt.MAXMOLTEMPLATEDIHEDRALS) then
          write(printstring,*) "Exceeded maximum number of molecule template dihedrals (",
     .                         MAXMOLTEMPLATEDIHEDRALS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXMOLTEMPLATEDIHEDRALS and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NMolTemplateDihedrals(itmpl) = idih
      endif
      
      MolTemplateDihName(itmpl,idih,1) = iatm1
      MolTemplateDihRes(itmpl,idih,1) = ires1
      MolTemplateDihName(itmpl,idih,2) = iatm2
      MolTemplateDihRes(itmpl,idih,2) = ires2
      MolTemplateDihName(itmpl,idih,3) = iatm3
      MolTemplateDihRes(itmpl,idih,3) = ires3
      MolTemplateDihName(itmpl,idih,4) = iatm4
      MolTemplateDihRes(itmpl,idih,4) = ires4
      MolTemplateDihedralType(itmpl,idih,1)= dummy
      MolTemplateDihedralType(itmpl,idih,2)= flex
      MolTemplateDihedralDelta(itmpl,idih) = deg2Rad(delta)
      MolTemplateDihPar(itmpl,idih,1) = ipar0
      MolTemplateDihPar(itmpl,idih,2) = ipar1
      MolTemplateDihNWells(itmpl,idih) = nwells
      MolTemplateDihWell(itmpl,idih,1) = well1
      MolTemplateDihWell(itmpl,idih,2) = well2
      MolTemplateDihWell(itmpl,idih,3) = well3
      MolTemplateDihWell(itmpl,idih,4) = well4
c      write(printstring,*) ' Tmpl ',itmpl,' dih ',idih
c      call printLine(DEBUG,printstring)
c      write(printstring,*) ' Nwells ',nwells,' well2 ',well2
c      call printLine(DEBUG,printstring)
      return
      end
