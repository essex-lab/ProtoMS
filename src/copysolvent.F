
      subroutine copySolvent(fromsvn,tosvn)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
c########################################################
c
c     Copy everything about solvent 'fromsvn' to 
c     solvent 'tosvn'. Note that both 'fromsvn' and
c     'tosvn' must be between 1 and NSolvents
c
c     (C) Christopher Woods, February 2005
c
c########################################################

      integer fromsvn,tosvn,i,j
      
      if (fromsvn.lt.1 .or. fromsvn.gt.NSolvents .or.
     .    tosvn.lt.1 .or. tosvn.gt.NSolvents) then
        write(printstring,*) "Trying to copy invalid solvents... ",fromsvn," ",tosvn
        call printLine(WARNING,printstring)
        return
      endif
      
c     copy the information
      NSvnAtoms(tosvn) = NSvnAtoms(fromsvn)
      
      do j=1,3
        do i=1,NSvnAtoms(fromsvn)
          ZSvnCoords(tosvn,i,j) = ZSvnCoords(fromsvn,i,j)
        enddo
        
        SvnCentre(tosvn,j) = SvnCentre(fromsvn,j)
      enddo
      
      SvnADel(tosvn) = SvnADel(fromsvn)
      SvnRDel(tosvn) = SvnRDel(fromsvn)
      
      SvnName(tosvn) = SvnName(fromsvn)
      SvnRes(tosvn) = SvnRes(fromsvn)
      
      do i=1,NSvnAtoms(fromsvn)
        ZSvnPar(tosvn,i) = ZSvnPar(fromsvn,i)
        SvnAtomName(tosvn,i) = SvnAtomName(fromsvn,i)
      enddo
      
      return
      end
