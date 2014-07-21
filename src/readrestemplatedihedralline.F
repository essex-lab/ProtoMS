      subroutine readResTemplateDihedralLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c########################################################
c
c     This routine reads in a residue dihedral template line
c
c     (C) Christopher Woods
c
c########################################################

      character*(*) line
      integer itmpl
      character*4 at1,at2,at3,at4
      logical flex,dummy
      integer ipar0,ipar1,idih,findTemplateDihedral
      double precision delta,deg2Rad
      
c     format is "dihedral atm1 atm2 atm3 atm4 (flex idelta) (dummy) (param ipar0 ipar1)

      if (nwords.lt.5) then
        call printLine(WARNING,"Incorrectly formatted residue dihedral template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)
      at3 = words(4)(1:4)
      at4 = words(5)(1:4)
            
      call sortTemplateDihedral(at1,at2,at3,at4)
      idih = findTemplateDihedral(at1,at2,at3,at4,itmpl)
      if (idih.eq.0) then
        idih = NTemplateDihedrals(itmpl) + 1
        if (idih.gt.MAXTEMPLATEDIHEDRALS) then
          write(printstring,*) "Exceeded maximum number of residue template dihedrals (",
     .                         MAXTEMPLATEDIHEDRALS,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXTEMPLATEDIHEDRALS and recompile")
          return
        endif 
        NTemplateDihedrals(itmpl) = idih
      else
        write(printstring,*) "Overwriting residue template angle for ",at1,"-",at2,
     .                      "-",at3,"-",at4
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
      
      TemplateDihName(itmpl,idih,1) = at1
      TemplateDihName(itmpl,idih,2) = at2
      TemplateDihName(itmpl,idih,3) = at3
      TemplateDihName(itmpl,idih,4) = at4
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      TemplateDihedralType(itmpl,idih,1) = dummy
      TemplateDihedralType(itmpl,idih,2) = flex
      TemplateDihedralDelta(itmpl,idih) = deg2Rad(delta)
      TemplateDihedralPar(itmpl,idih,1) = ipar0
      TemplateDihedralPar(itmpl,idih,2) = ipar1
      
      return
      end
