      subroutine getDihWells(well1,well2,well3,well4,nwells)
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'printer.inc'
c##################################################################
c
c     This routine scans the words to find out whether a 
c     a dihedral has several defined 'wells'
c
c     Julien Michel March 05
c##################################################################

      integer i
      double precision getDouble
      integer nwells
      double precision well1,well2,well3,well4
      double precision deg2rad
      well1 = 0.0d0
      well2 = 0.0d0
      well3 = 0.0d0
      well4 = 0.0d0
      nwells = 0
c     scan through the words to find 'dummy', 'flex' or 'param'      
      do i=1,nwords
        if (words(i) .eq. 'wells') then
          call split2(words(i+1),",")
          well1 = getDouble(words2(1),0.0)
          well1 = deg2rad(well1)
          well2 = getDouble(words2(2),0.0)
          well2 = deg2rad(well2)
          well3 = getDouble(words2(3),0.0)
          well3 = deg2rad(well3)
          well4 = getDouble(words2(4),0.0)
          well4 = deg2rad(well4)
          nwells = nwords2
c          write(printstring,*) 'Wells for this Dihedral ',well1,well2,well3,well4,nwells
c          call printLine(DEBUG,printstring)
          return
        endif
      enddo
      
      return
      end
