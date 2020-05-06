      subroutine readResTemplateBondLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c#######################################################
c
c     This routine reads in residue template bonds
c
c     (C) Christopher Woods  September 2004
c
c#######################################################
            
      character*(*) line
      integer itmpl,ibnd,findTemplateBond
      character*4 at1,at2
      logical flex,dummy
      integer ipar0,ipar1
      double precision delta
      
c     format is "bond atm1 atm2 (flex idelta) (dummy) (param ipar0 ipar1)"
c     (items in parenthesis are optional)
      
      if (nwords.lt.3) then
        call printLine(WARNING,"Incorrectly formatted residue bond template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)

      call sortTemplateBond(at1,at2)
      
      ibnd = findTemplateBond(at1,at2,itmpl)
      if (ibnd.eq.0) then
        ibnd = NTemplateBonds(itmpl) + 1
        if (ibnd.gt.MAXTEMPLATEBONDS) then
           write(printstring,*) "Exceeded maximum number of residue template bonds (",
     .                           MAXTEMPLATEBONDS,")"
           call printLine(FATAL,printstring)
           write(printstring,*)"Increase MAXTEMPLATEBONDS and recompile"
           call closeProgram(-1,printstring)
           return
        endif
        NTemplateBonds(itmpl) = ibnd
      else                          
        write(printstring,*) "Overwriting residue template bond for ",at1,"-",at2
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
      
      TemplateBndName(itmpl,ibnd,1) = at1
      TemplateBndName(itmpl,ibnd,2) = at2
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      TemplateBondType(itmpl,ibnd,1) = dummy
      TemplateBondType(itmpl,ibnd,2) = flex
      TemplateBondDelta(itmpl,ibnd) = delta
      TemplateBondPar(itmpl,ibnd,1) = ipar0
      TemplateBondPar(itmpl,ibnd,2) = ipar1
      
      return      

      end  
