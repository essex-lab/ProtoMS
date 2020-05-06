
c  This file provides the variables used in the move lists

c     the maximum number of residues that have moved that can be stored
c     in the residue move list
      integer MAXRESMOVED
      parameter(MAXRESMOVED=500)
      
c     the maximum number of solutes that have moved that can be stored
      integer MAXSOLMOVED
      parameter(MAXSOLMOVED=20)
      
c     the maximum number of solvents that have moved that can be stored
      integer MAXSVNMOVED
      parameter(MAXSVNMOVED=50)
      
c     list of all moved residues, (protein,residue) numbers
      integer resMovedList(MAXRESMOVED,2)
      
c     list of all moved solutes
      integer solMovedList(MAXSOLMOVED)
      
c     list of all moved solvents
      integer svnMovedList(MAXSVNMOVED)
      integer gcsMovedList(MAXSVNMOVED)
c     whether or not everything has moved
      logical allmoved
      
c     the number of moved residues,solutes,solvents
      integer nResMoved,nSolMoved,nSvnMoved,nGCSMoved
      
      common /movelists/ resMovedList,solMovedList,svnMovedList,nResMoved,
     .                   nSolMoved,nSvnMoved,allmoved,gcsMovedList,nGCSMoved
     
      save /movelists/
