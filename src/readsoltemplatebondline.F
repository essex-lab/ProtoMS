      subroutine readSolTemplateBondLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###########################################################
c
c     This routine reads in the bonds for a molecule template
c
c     (C) Christopher Woods
c
c###########################################################

      character*(*) line
      integer itmpl
      character*4 iatm1,ires1,iatm2,ires2
      logical dummy,flex
      double precision delta
      integer ipar0,ipar1
      integer findMolTemplateBond
      integer ibnd
      
c     format is 'bond iatm1 ires1 iatm2 ires2 (dummy) (flex delta) (param ipar0 ipar1)'
c      (items in parenthesis are optional)
      
      if (nwords.lt.5) then
        call printLine(WARNING,"Incorrectly formatted solute template bond line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      iatm1 = words(2)(1:4)
      ires1 = words(3)(1:4)
      iatm2 = words(4)(1:4)
      ires2 = words(5)(1:4)
      
c     get the optional parameters
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
        
c     see if the bond already exists
      ibnd = findMolTemplateBond(iatm1,ires1,iatm2,ires2,itmpl)      
      if (ibnd.ne.0) then
        write(printstring,*) "Overwriting information for solute bond (",
     .         iatm1,",",ires1,")-(",iatm2,",",ires2,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      else
        ibnd = NMolTemplateBonds(itmpl) + 1
        if (ibnd.gt.MAXMOLTEMPLATEBONDS) then
          write(printstring,*) "Exceeded maximum number of molecule template bonds (",
     .                         MAXMOLTEMPLATEBONDS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXMOLTEMPLATEBONDS and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NMolTemplateBonds(itmpl) = ibnd
      endif
      
      MolTemplateBndName(itmpl,ibnd,1) = iatm1
      MolTemplateBndRes(itmpl,ibnd,1) = ires1
      MolTemplateBndName(itmpl,ibnd,2) = iatm2
      MolTemplateBndRes(itmpl,ibnd,2) = ires2
      MolTemplateBondType(itmpl,ibnd,1)= dummy
      MolTemplateBondType(itmpl,ibnd,2)= flex
      MolTemplateBondDelta(itmpl,ibnd) = delta
      MolTemplateBndPar(itmpl,ibnd,1) = ipar0
      MolTemplateBndPar(itmpl,ibnd,2) = ipar1
      
      return
      end
