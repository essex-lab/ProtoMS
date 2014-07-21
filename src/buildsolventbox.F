
      subroutine buildSolventBox(xmin,ymin,zmin,xmax,ymax,zmax)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
c###################################################################
c
c     Build a solvent box with from (xmin,ymin,zmin) to 
c     (xmax,ymax,zmax). 
c
c     (C) Christopher Woods, February 2005
c
c###################################################################

      double precision xmin,ymin,zmin,xmax,ymax,zmax

      double precision vec(3),smin(3),smax(3),sdim(3),rmin(3),rmax(3),dist2
      double precision icoor(3),irad2
      integer imin(3),imax(3)
      integer i,j,k,isvn,nsvn
      integer iprot,isol,iat
      
      if (NSolvents.le.0) then
        call printLine(WARNING,"Cannot replicate a solvent box if no solvents are loaded!")
        return
      endif
      
      rmin(1) = xmin
      rmin(2) = ymin
      rmin(3) = zmin
      rmax(1) = xmax
      rmax(2) = ymax
      rmax(3) = zmax

      write(printstring,*) (rmin(i),i=1,3),(rmax(i),i=1,3)
      call printLine(INFO,printstring)
      
c     This is a very quick and dirty algorithm! We assume that the currently
c     loaded solvent molecules are in a dense cuboid, which we will then 
c     replicate to fill the whole volume. We will then remove the solvent molecules
c     that are too close to any solutes. This is a very wasteful algorithm as
c     it will add too many solvent molecules (into the volume occupied by the
c     solute) and is just not good! I really need to have a better C/C++ function
c     to do this...

c     get the dimensions of the current box from whatever the user has set
      
c     make sure something has been set, otherwise calculate some...
      if (BoxDimension(1).eq.ZERO) then
        call calculateSolventBox()
      endif
      
      do i=1,3
        smin(i) = BoxOrigin(i)
        smax(i) = BoxTop(i)
        sdim(i) = BoxDimension(i)
      enddo
      
      write(printstring,*) "Existing solvent box extends from (",smin(1),",",smin(2),",",smin(3),") to (",
     .        smax(1),",",smax(2),",",smax(3),"), dimension (",sdim(1),",",sdim(2),",",sdim(3),")"
      call printLine(INFO,printstring)
      
c     now work out how many replications we will need
      do i=1,3
        if (rmin(i).ge.smin(i)) then
          imin(i) = 0
        else
          imin(i) = int( (rmin(i) - smin(i)) / sdim(i) ) - 1
        endif
        
        if (rmax(i).le.smax(i)) then
          imax(i) = 0
        else
          imax(i) = int( (rmax(i) - smax(i)) / sdim(i) ) + 1
        endif
      enddo
      
c     now loop over all of the replicas
      write(printstring,*) (imin(i),i=1,3),(imax(i),i=1,3)
      call printLine(INFO,printstring)
      
      call printLine(INFO,"Replicating solvent box...")
      
c     save the number of solvents before replication
      nsvn = NSolvents      
      
      do i=imin(1),imax(1)
        vec(1) = i*sdim(1)
        
        do j=imin(2),imax(2)
          vec(2) = j*sdim(2)
          
          do k=imin(3),imax(3)
            if (.not. (i.eq.0 .and. j.eq.0 .and. k.eq.0) ) then
              vec(3) = k*sdim(3)
              
c             loop over all solvent molecules...
              do isvn=1,nsvn
c               create a copy of this solvent molecule
                NSolvents = NSolvents + 1
                if (NSolvents.gt.MAXSOLVENTS) then
                  call closeProgram(-1,"Need to increase MAXSOLVENTS and recompile!")
                  return
                endif
                 
                call copySolvent(isvn,NSolvents)
c               translate the solvent to its new location
                call translateSolvent(NSolvents,vec)
              enddo
            endif
          enddo
        enddo
      enddo                                        

c     now remove any solvent molecules that are outside of the box                                
      nsvn = NSolvents
      
      do i=nsvn,1,-1
        if (SvnCentre(i,1).gt.xmax .or. SvnCentre(i,1).lt.xmin .or.
     .      SvnCentre(i,2).gt.ymax .or. SvnCentre(i,2).lt.ymin .or.
     .      SvnCentre(i,3).gt.zmax .or. SvnCentre(i,3).lt.zmin) then
c         delete the solvent molecule            
          if (i.eq.NSolvents) then
            NSolvents = NSolvents - 1
          else
c           copy the last solvent to this position
            call copySolvent(NSolvents,i)
            NSolvents = NSolvents - 1
          endif
        endif
      enddo            

      nsvn = NSolvents
      
c     remove any close contacts...
      do i=nsvn,1,-1

c       copy the center of the solvent into a tmp array      
        do k=1,3
          icoor(k) = SvnCentre(i,k)
        enddo
      
c       get the solvent radius squared
        irad2 = ZERO
        do iat=1,NSvnAtoms(i)
          dist2 = (ZSvnCoords(i,iat,1)-icoor(1))**2 + (ZSvnCoords(i,iat,2)-icoor(2))**2 + 
     .                      (ZSvnCoords(i,iat,3)-icoor(3))**2
          
          if (dist2.gt.irad2) irad2 = dist2
        enddo
        
c       add a couple of A to the solvent radius to account for the volume of the atoms
        irad2 = irad2 + FOUR        
        
c       loop over all protein atoms and check for close contacts
        do iprot=1,NProteins
          do iat=1,NZAtoms(iprot)
            dist2 = (ZCoords(iprot,iat,1)-icoor(1))**2 + (ZCoords(iprot,iat,2)-icoor(2))**2 + 
     .                      (ZCoords(iprot,iat,3)-icoor(3))**2
            
c           if this distance is less than the cutoff then delete this solvent molecule            
            if (dist2.lt.irad2) goto 100
          enddo
        enddo
        
c       loop over all solute atoms and check for close contacts
        do isol=1,NSolutes
          do iat=1,NSolAtoms(isol)
            dist2 = (ZSolCoords(isol,iat,1)-icoor(1))**2 + (ZSolCoords(isol,iat,2)-icoor(2))**2 +
     .                  (ZSolCoords(isol,iat,3)-icoor(3))**2
                   
c           if this distance is less than the cutoff then delete this solvent molecule            
            if (dist2.lt.irad2) goto 100
          enddo
        enddo

        goto 110
            
100     continue
c       delete this solvent molecule        
        if (i.eq.NSolvents) then
          NSolvents = NSolvents - 1
        else
          call copySolvent(NSolvents,i)
          NSolvents = NSolvents - 1
        endif

110     continue
c       don't delete the solvent molecule
      enddo

c     now set the new solvent box dimensions
      do i=1,3
        BoxOrigin(i) = rmin(i)
        BoxTop(i) = rmax(i)
        BoxDimension(i) = rmax(i) - rmin(i)
        BoxDimension2(i) = HALF*BoxDimension(i)
      enddo
            
      return
      end
      
