
      subroutine processSolventLine(check)
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'printer.inc'
c###################################################
c
c     This routine processes a 'solvent' chunk. 
c     These chunks are used to build solvent boxes
c     or solvent caps.
c
c     Note that these chunks only do basic work
c     of building the boxes/caps. To get a high
c     quality box/cap you should perform NPT
c     equilibration of the box or NVT equilibration
c     of the cap with included solvent swap moves.
c
c     If 'check' is true, then only print out
c     what we will do, don't actually do it.
c
c     (C) Christopher Woods, February 2005
c
c###################################################

      logical check
      double precision xmin,ymin,zmin,xmax,ymax,zmax,itmp
      double precision ox,oy,oz,irad
      double precision getDouble
      
c     The following are commands used to build a solvent box/cap. The box
c     is made by replicating the currently loaded solvent until it reaches
c     the correct dimensions (or removing solvent if the box/cap is smaller
c     than that already loaded). 
c
c     This is then followed by an elimination phase where solvents that are
c     too close to another molecule are deleted.
c
c     Finally, an addition phase scans the box/cap for any holes sufficiently
c     large to hold a solvent, and adds one in there at a random orientation.
c     This is likely to be quite an expensive operation!
c
c     #build a solvent box centered on 0,0,0 with dimensions (xdim,ydim,zdim)
c     solvent box xdim ydim zdim                   
c     #build a solvent box with origin (xorig,yorig,zorig) going to (xmax,ymax,zmax)
c     solvent box xorig yorig zorig xmax ymax zmax
c
c     #build a solvent cap centered at (xorig,yorig,zorig) with radius 'rad'
c     solvent cap xorig yorig zorig rad

      if (nwords.lt.5) then
        write(printstring,*) "Incorrent format for solvent chunk line"
        call printLine(WARNING,printstring)
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return
      endif
      
c     see whether we are building a box or a cap
      call lowercase(words(2))      
      
      if (words(2).eq.'box') then
c       we are building a box - now read the box dimensions
        if (nwords.lt.8) then
c         read in xdim,ydim,zdim
          xmin = getDouble(words(3),ZERO)
          ymin = getDouble(words(4),ZERO)
          zmin = getDouble(words(5),ZERO)
          
          write(printstring,*) "Creating a solvent box centred on the origin with dimensions (",
     .             xmin,",",ymin,",",zmin,")"
          call printLine(INFO,printstring)
          
c         convert to xmin/xmax
          xmax = HALF * dabs(xmin)
          xmin = -xmax
          ymax = HALF * dabs(ymin)
          ymin = -ymax
          zmax = HALF * dabs(zmin)
          zmin = -zmax
        else
c         read the complete set of dimensions
          xmin = getDouble(words(3),ZERO)
          ymin = getDouble(words(4),ZERO)
          zmin = getDouble(words(5),ZERO)
          xmax = getDouble(words(6),ZERO)
          ymax = getDouble(words(7),ZERO)
          zmax = getDouble(words(8),ZERO)
          
c         make sure that the dimensions are sane
          if (xmin .gt. xmax) then
            itmp = xmin
            xmin = xmax
            xmax = itmp
          endif
          
          if (ymin .gt. ymax) then
            itmp = ymin
            ymin = ymax
            ymax = itmp
          endif
          
          if (zmin .gt. zmax) then
            itmp = zmin
            zmin = zmax
            zmax = itmp
          endif

          write(printstring,*) "Building a solvent box from (",xmin,",",ymin,",",zmin,") to (",
     .                                    xmax,",",ymax,",",zmax,")"
          call printLine(INFO,printstring)
        endif

c       now build the solvent box...
        if (.not.check) call buildSolventBox(xmin,ymin,zmin,xmax,ymax,zmax)
                  
      else if (words(2).eq.'cap') then
c       we are building a solvent cap - read in the cap dimensions
        if (nwords.lt.6) then
          call printLine(WARNING,"Cannot interpret build solvent cap line!")
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif

c       read the origin of the solvent cap sphere                
        ox = getDouble(words(3),ZERO)
        oy = getDouble(words(4),ZERO)
        oz = getDouble(words(5),ZERO)

c       read the radius of the sphere
        irad = dabs(getDouble(words(6),ZERO))
        
c       get the box dimensions that hold this sphere
        xmin = ox - irad
        xmax = ox + irad
        ymin = oy - irad
        ymax = oy + irad
        zmin = oz - irad
        zmax = oz + irad        
        
        write(printstring,*) "Building a solvent cap centered at (",ox,",",oy,",",oz,"), radius ",
     .                       irad
        call printLine(INFO,printstring)
        
c       build a box that contains the cap
        if (.not.check) then
          call buildSolventBox(xmin,ymin,zmin,xmax,ymax,zmax)

c         now remove any solvent molecules that lie outside the cap
          call removeExcessCapSolvent(ox,oy,oz,irad)
        endif

      endif      
      
      return
      
      end
      
