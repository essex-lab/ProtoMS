      subroutine syncSoluteDummies(isol1,isol2)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'enums.inc'
c#############################################################
c
c     This routine sets solutes isol1 and isol2 to have the
c     same dummy atoms. This is useful if you want to ensure that
c     the two solutes are translated and rotated about the same
c     points (as solute translation and rotation are about their
c     dummy atoms).
c
c     (C) Christopher Woods
c 
c#############################################################

      double precision ivec(3)
      integer isol1,isol2,i,j

      call prettyLine(INFO,"Syncing the dummy atoms of a pair of solutes")

      write(printstring,10) isol1,(SolDummy(isol1,1,i),i=1,3),
     .   (SolDummy(isol1,2,i),i=1,3),(SolDummy(isol1,3,i),i=1,3)
      call printLine(INFO,printstring)
      write(printstring,10) isol2,(SolDummy(isol2,1,i),i=1,3),
     .   (SolDummy(isol2,2,i),i=1,3),(SolDummy(isol2,3,i),i=1,3)
      call printLine(INFO,printstring)

10    format('Solute ',i3,' dummies =',3(' (',f8.3,',',f8.3,',',f8.3,')'))

c     move the dummies so that they are equidistant between the
c     dummies of both solutes (for want of something better!)
      do i=1,3
        ivec(i) = HALF * (SolDummy(isol2,1,i) - SolDummy(isol1,1,i))
      
        do j=1,3
          SolDummy(isol1,j,i) = SolDummy(isol1,j,i) + ivec(i)
          SolDummy(isol2,j,i) = SolDummy(isol1,j,i)
        enddo
      enddo

c     as the dummies have changed, we need to rebuild the zmatrix
      call getSoluteZmatrix(isol1)
      call getSoluteZmatrix(isol2)
      
      write(printstring,11) (SolDummy(isol2,1,i),i=1,3),
     .   (SolDummy(isol2,2,i),i=1,3),(SolDummy(isol2,3,i),i=1,3)
      call printLine(INFO,printstring)
      
11    format('New dummies =',3(' (',f8.3,',',f8.3,',',f8.3,')'))
      
      end
