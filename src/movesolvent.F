      subroutine moveSolvent(isvn,ienew)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
c###################################################
c
c     This routine moves solvent molecule isvn
c     via translation and rotation
c
c     Translation and rotation is done in the 
c     solvents local vector space
c
c     (C) Christopher Woods   4 March 2002
c     Julien Michel October 2004
c     MSB 2012 - prevents solvent from entering a GC
c     delimited region
c     New way to perform translations/rotations
c###################################################

      integer isvn,i,j,nats
      double precision ranRange,ienew
      double precision ranval,ranvec(3),idx(3)
      double precision distance

c     initialise the value of ienew to 0
      ienew = 0

c      PRINT*,"old" ,SvnCentre(isvn,1) ,SvnCentre(isvn,2) ,SvnCentre(isvn,3)

c     translate the solvent molecule randomly
      call getVectorOnSphere(ranvec)
      ranval = ranRange(-SvnADel(isvn),SvnADel(isvn))
      do i=1,3
        ranvec(i) = ranval * ranvec(i)
      enddo
      call translateSolvent(isvn,ranVec)
c     now rotate the solvent
      call getVectorOnSphere(ranvec)
      ranval = ranRange(-SvnRDel(isvn),SvnRDel(isvn))
      call rotateSolvent(isvn,ranvec,ranval)
c      PRINT*,"Post moved" ,SvnCentre(isvn,1) ,SvnCentre(isvn,2) ,SvnCentre(isvn,3)

c     need to map back into periodic box (if necessary)
      if (boundary.eq.PERIODICBOUNDARIES) then
        nats = NSvnAtoms(isvn)
c       if the center of the solvent is outside the box, then map it
c       over to the other side
        do i=1,3
          if (SvnCentre(isvn,i).gt.BoxTop(i)) then
            idx(i) = -BoxDimension(i)
          else if (SvnCentre(isvn,i).lt.BoxOrigin(i)) then
            idx(i) = BoxDimension(i)
          else 
            idx(i) = ZERO
          endif
          if (idx(i).ne.ZERO) then
            do j=1,nats
              ZSvnCoords(isvn,j,i) = ZSvnCoords(isvn,j,i) + idx(i)
            enddo
            SvnCentre(isvn,i) = SvnCentre(isvn,i) + idx(i)
          endif
        enddo
      endif


      if (sphericalgcmc) then
      distance = ((SvnCentre(isvn,1)-centerx)**2.+(SvnCentre(isvn,2)
     :  -centery)**2+(SvnCentre(isvn,3)-centerz)**2)**0.5

      if (distance.lt.gcradius) then
      print*, 'movesolvent'
      goto 99
      endif
      endif

      if (SvnCentre(isvn,1).gt.(coordx+gridx)) then
      goto 99
      else
      goto 1
      endif
1      if (SvnCentre(isvn,2).gt.(coordy+gridy)) then
      goto 99
       else
       goto 2
      endif
2      if (SvnCentre(isvn,3).gt.(coordz+gridz)) then
      goto 99
       else
       goto 3
      endif
3      if (SvnCentre(isvn,1).lt.(coordx)) then
      goto 99
      else
      goto 4
      endif
4      if (SvnCentre(isvn,2).lt.(coordy)) then
      goto 99
      else
      goto 5
      endif
5      if (SvnCentre(isvn,3).lt.(coordz)) then
       goto 99
       else
      goto 6
      endif

6      ienew = (ienew + (HUGE*HUGE))
      goto 99
99    continue
      return
      end
