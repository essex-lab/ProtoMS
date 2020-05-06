      double precision function getZmatDihedral(ityp,imol,ires,iatm,dihf,dihb)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
c################################################
c     Return the value of a a dihedral angle that is
c     specified in a zmatrix
c
c     Julien Michel February 2005
c
c#################################################
      integer ityp,imol,ires,iatm
      double precision dihf,dihb
      
      getZmatDihedral = 0.0d0
      if ( ityp .eq. SOLUTE) then
        getZmatDihedral = ZSolVal(imol,iatm,3)
        if (NSolVariable(imol) .gt. 0) then 
          dihf = ZSolVal(NSolForwards(imol),iatm,3)
          dihb = ZSolVal(NSolBackwards(imol),iatm,3)
        else
          dihf = getZmatDihedral
          dihb = getZmatDihedral
        endif
      else if ( ityp .eq. PROTEIN) then
        getZmatDihedral = ZValSC(imol,ires,iatm,3)
        dihf = getZmatDihedral
        dihb = getZmatDihedral
      else 
        write(printstring,*) 'Requested Zmat torsion of ',
     .    ' an atom that does not have a zmatrix !'
        call printLine(FATAL,printstring)
        call closeProgram(-1)
      endif

      return
      end

       
