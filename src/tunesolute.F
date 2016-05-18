      double precision function tuneSolute(isol,ires,idih)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'printer.inc'
c#######################################################
c
c     This routine moves a single residue (ires) of
c     solute isol.
c
c     (C) Christopher Woods    28 February 2002
c     Rewritten by Julien Michel October 2004 
c#######################################################
      
      integer isol,ires,i,j
      integer istrt,iend
      double precision ranRange
      double precision ranval,ranvec(3)
      integer ranInt
      integer jsol,jres
      integer ranwell
      integer iDih,jDih
c$$$      print *, idih
c
c     Determine if this is a dual topology solute
c
      jsol = DualTopologySolute(isol)
      if (jsol .lt. 0) jsol = -jsol

c     In group dual topology, each solute is moved independently
      if (NSolGroups.gt.ZERO) jsol = 0

c
c     Translation 
c
      if (SolADel(isol).gt.ZERO) then 
        call getVectorOnSphere(ranvec)
        ranval = ranRange(-SolADel(isol),SolADel(isol))
        do i=1,3
          ranvec(i) = ranval * ranvec(i)
        enddo
        call translateSolute(isol,ranVec)
      endif

      if (jsol .ne. 0) then
        if (DualTopologySyncTrans(isol)) then
          if (SolADel(jsol) .gt. ZERO .and. SolADel(isol).gt.ZERO) call translateSolute(jsol,ranvec)
        else
          call getVectorOnSphere(ranvec)
          ranval = ranRange(-SolADel(jsol),SolADel(jsol))
          do i=1,3
            ranvec(i) = ranval * ranvec(i)
          enddo
          call translateSolute(jsol,ranvec)
        endif
      endif

c
c     Rotation
c
      if (SolRDel(isol).gt.ZERO) then
        call getVectorOnSphere(ranvec)
        ranval = ranRange(-SolRDel(isol),SolRDel(isol))
        call rotateSolute(isol,ranvec,ranval)
      endif

      if (jsol .ne. 0) then
        if (DualTopologySyncRot(isol)) then
          if (SolRdel(jsol) .gt. ZERO .and. SolRdel(isol).gt.ZERO) call rotateSolute(jsol,ranvec,ranval)
        else
          call getVectorOnSphere(ranvec)
          ranval = ranRange(-SolRDel(jsol),SolRDel(jsol))
          call rotateSolute(jsol,ranvec,ranval)
        endif
      endif

c
c     Internal Coordinates
c

c     find the range of atoms encompased by the residue

      istrt = SolResidues(isol,ires,1)
      iend = SolResidues(isol,ires,2)
      
      jDih = 0
c     loop over the range of atoms and move them
      do i=istrt,iend
        do j=1,3
           if (ZSolMov(isol,i,j)) then
              if (j.eq.3) then
                 jDih = jDih + 1
c$$$                 print *, ZSolVal(isol,i,j)
                 if (jDih.eq.iDih) then

                    if (ZSolDihNWells(isol,i).gt.0) then
                       ranwell = ranint(1,ZSolDihNWells(isol,i))
                       ranval = ranval + ZSolDihWell(isol,i,ranwell)
                    else
                       ranval = ranRange(-ZSolDelta(isol,i,j),ZSolDelta(isol,i,j))
                    endif
                    tunesolute = ranval
                 else
                    ranval = 0
                 endif
c$$$                 print *,'moving', i, ranval
              else
                 ranval = ranRange(-ZSolDelta(isol,i,j),ZSolDelta(isol,i,j))
              endif
              ZSolVal(isol,i,j) = ZSolVal(isol,i,j) + ranval
           endif
         enddo
      enddo

      if (jsol .ne. 0) then 
c       We need to pick one residue...
        jres = ranInt(1,NSolResidues(jsol))
        istrt = SolResidues(jsol,jres,1)
        iend = SolResidues(jsol,jres,2)
        do i=istrt,iend
          do j=1,3
            if (ZSolMov(jsol,i,j)) then
              ranval = ranRange(-ZSolDelta(jsol,i,j),ZSolDelta(jsol,i,j))
c           If there are wells defined for this dihedral, randomly pick
c           one well
              if (j .eq. 3 .and. ZSolDihNWells(isol,i) .gt. 0) then
                ranwell = ranint(1,ZSolDihNWells(isol,i))
                ranval = ranval + ZSolDihWell(isol,i,ranwell)
              endif
              ZSolVal(jsol,i,j) = ZSolVal(jsol,i,j) + ranval
            endif
          enddo
        enddo        
      endif

      return
      end
