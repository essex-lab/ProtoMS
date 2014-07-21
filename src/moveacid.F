      subroutine moveAcid(iprot,ires)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'printer.inc'
C     This is the funky subroutine that moves the amino
C     acid backbone, one residue at a time. Residue
C     ires in protein iprot is moved, and 
C     movetype and rreqdel can be used to 
C     specify the move type and change in move...

c     changed to move the N-C-CA atoms together, as
c     this should better preserver the geometry of the
c     aminoacids...

c     #################
c     # Historical note - this was
c     # one of the first routines written
c     # for ProtoMS, and this move
c     # started it all...
c     #################

C     movetype          action
C     0                 Random moves + random deltas
C     1,2,3             Translation in x,y or z, by rreqdel
C     4                 Pincer move by rreqdel
C     5,6,7             Rotation around x,y or z, by rreqdel

C     (C) Christopher Woods  16 January 2002

      double precision vecx(3),vecCO(3),vecy(3),vecz(3)
      double precision ori(3),vecCaC(3),vecNC(3)
      double precision veclength2
      double precision beta
      double precision dCaN,dCaN2,dCaC2,dNC2
      double precision rotdel
      double precision trandel
      double precision del
      integer iprot,ires,ires1,ires2
      integer i,j
      double precision ranval,ranRange
c#########################################

c     moving N2 C1 CA1 as a unit (which will also move HN2,O1,HA1,CB1)

c       O2  HA2 HN2      CB1
c       ||   |   |        |
c     --C2--CA2--N2--C1--CA1--N1--
c            |       ||   |   |
c                    O1  HA1  HN1

c     since the move involves atoms from the next residue, we must
c     not be on the last residue
      if (ires.gt.NResidues(iprot).or.ires.le.0) then
        write(printstring,*) "Trying to move an invalid residue ",ires,NResidues(iprot)
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      else if (ires.eq.NResidues(iprot)) then
c       this is the last residue - since the move is with this and the
c       next residue, we need to do something different - we instead 
c       translate the last three backbone atoms randomly (since the first
c       atom will be moved by the previous residue)
        trandel = ResADel(iprot,ires)*0.1
        
        do i=1,3
c         translate the backbone          
          do j=2,4
c           pick a random amount to translate by
            ranval = ranRange(-trandel,trandel)
            BBCoord(iprot,ires,j,i) = BBCoord(iprot,ires,j,i) + ranval
          enddo
        enddo
        do i=1,3
          do j=1,4
            ZCoords(iprot,BBZid(iprot,ires,j),i) = BBCoord(iprot,ires,j,i)
          enddo
        enddo
        return
      endif 
            
      ires1 = ires
      ires2 = ires + 1

c     take the current positions of the system, and
c     get the vector between CA1 and N2... (+ve towards N2)
      do i=1,3
        vecx(i) = BBCoord(iprot,ires2,1,i) - BBCoord(iprot,ires1,2,i)
      enddo
c     get the length of the CA-N bond (will need it later)
      dCaN2 = veclength2(vecx)
      dCaN = sqrt(dCaN2)
c     normalise the vector
      call vecnormalise(vecx)

c     now find the vector of the carbonyl bond... (+ve towards O1)
      do i=1,3
        vecCO(i) = BBCoord(iprot,ires1,4,i)-BBCoord(iprot,ires1,3,i)
      enddo
      call vecnormalise(vecCO)

c     get the x/y/z vectors based on vecx and vecCo      
      call getAxisVectors(vecx,vecCO,vecy,vecz)
      
c     now that we the x/y/z vectors,
c     it is time to start moving the molecule. First we translate
c     the whole amino acid by random amounts along vecx,vecy,vecz
c     translating the molecule - random translation about all three
c     axes if movetype = 0, else ordered translation about one of 
c     the axis if requested
      trandel = ResADel(iprot,ires)
c     get a random dble from -xtrandel to xtrandel
      del = ranRange(-trandel,trandel)
      call translateAcid(iprot,ires,vecx,del)
      del = ranRange(-trandel,trandel)
      call translateAcid(iprot,ires,vecy,del)
      del = ranRange(-trandel,trandel)
      call translateAcid(iprot,ires,vecz,del)

c     if this is the first residue, then the terminal N has not been
c     moved - translate it randomly (eventually should build and move
c     nitrogen via implicit zmatrix via BBAtoms 2-4)
      if (ires1.eq.1) then
        trandel = ResADel(iprot,ires) * 0.25
        do i=1,3
          del = ranRange(-trandel,trandel)
          BBCoord(iprot,ires1,1,i) = BBCoord(iprot,ires1,1,i) + del
        enddo
      endif        

c     We must now find the origin of rotation - this is 
c     the point along vecx that is directly (along vecy) 
c     below the CO bond
      do i=1,3
        vecCaC(i) = BBCoord(iprot,ires1,2,i) - BBCoord(iprot,ires1,3,i)
        vecNC(i) = BBCoord(iprot,ires2,1,i) - BBCoord(iprot,ires1,3,i)
      enddo
      dCaC2 = veclength2(vecCaC)
      dNC2 = veclength2(vecNC)
            
      beta = (dCaC2-dNC2-dCaN2) / (2*dCaN)
      do i=1,3
        ori(i) = BBCoord(iprot,ires2,1,i) + (beta*vecx(i))
      enddo

c     now we have the origin of rotation, rotate around the three
c     principle axes
      rotdel = ResRDel(iprot,ires)
c     now we can rotate the acid around vecx and ori by xdel radians
c     generate a random dble
      del = ranRange(-rotdel,rotdel)
      call rotateAcid(iprot,ires,ori,vecx,del)
c     do the same for vecy and vecz
      del = ranRange(-rotdel,rotdel)
      call rotateAcid(iprot,ires,ori,vecy,del)
      del = ranRange(-rotdel,rotdel)
      call rotateAcid(iprot,ires,ori,vecz,del)

c     now that the backbone has been moved, synch the coordinates
c     from residue to zmatrix based arrays (C,CA,N) with the 
c     zmatrix-id indexed array (ZCoords)
      do i=1,3
        do j=1,4
          ZCoords(iprot,BBZid(iprot,ires1,j),i) = BBCoord(iprot,ires1,j,i)
          ZCoords(iprot,BBZid(iprot,ires2,j),i) = BBCoord(iprot,ires2,j,i)
        enddo
      enddo

      return
      end
