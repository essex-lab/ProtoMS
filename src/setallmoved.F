      subroutine setAllMoved
      include 'dimensions.inc'
      include 'movelist.inc'
c###########################################
c
c     This is used to set everything as
c     having moved
c
c     (C) Christopher Woods, October 2004
c
c###########################################

      allMoved = .true.
      
c     clear the other lists as they are now unnecessary
      nResMoved = 0
      nSolMoved = 0
      nSvnMoved = 0
      nGCSMoved = 0
      return
      end
