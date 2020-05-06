      subroutine clearid()
      include 'dimensions.inc'
      include 'labels.inc'
      include 'splitter.inc'
      include 'printer.inc'     
c#############################################
c     
c     This routine clear an assigned ID number
c     
c     The format of an clearid line is :
c
c     chunkXX id clear 1 
c     
c     Julien Michel, February 2005
c
c############################################
      integer id
      integer getInt

      if (nwords .ne. 3) then
        write(printstring,*) 'Could not interpret clearid line !'
        call printLine(WARNING,printstring)
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
      endif
      id = getInt(words(3),0)
      write(printstring,*) 'Clearing ID number ',id
      call printLine(INFO,printstring)
      
      if (LabelID(id) .eq. 0) then
        write(printstring,*) 'ID ',id,' was not set ?'
        call printLine(WARNING,printstring)
      endif

      LabelID(id) = 0
      LabelType(id) = 0
      LabelMolecule(id) = 0
      LabelResidueName(id) = NULLNAME
      LabelResidueNumber(id) = 0
      LabelName(id) = NULLNAME

      return 
      end
