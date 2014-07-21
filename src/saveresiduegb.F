      subroutine SaveResidueGB(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'printer.inc'
      include 'gb.inc'
c############################################
c     Save residue coordinates on to a list
c     of moving residues for GB calculations
c
c     Julien Michel November 2004
c##############################################
      integer iprot,ires,i,j,k,resstart,resend

      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return
      
      countGBResidue = countGBResidue + 1
      if (countGBResidue .gt. MAXRESIDUES) then
        call printLine(WARNING,"Increase MAXRESIDUES  and recompile")
        call closeProgram(-1,"ProtoMC exiting !")
      endif
      resstart = BBZid(iprot,ires,1)
      resend = BBZid(iprot,ires,4)

c     We must copy the way coordinates in Zcoords are saved, but we
c     want to start from 1, not resstart

c     now save the coordinates
      k = 1
      do i=resstart,resend
        gbResidueIndex(iprot,i) = k
        do j=1,3
          gbResidueCoords(countGBResidue,k,j) = ZCoords(iprot,i,j)
        enddo
        k = k + 1
      enddo
c     DEBUG
c      k = 1
c      do i=resstart,resend
c        write(printstring,*) ' i ',i,' Coords ',ZCoords(iprot,i,1),' ',
c     &    ZCoords(iprot,i,2),' ',ZCoords(iprot,i,3)
c        call printLine(DEBUG,printstring)
c        write(printstring,*) ' k ',k,' gb coords ',
c     &    gbResidueCoords(countGBResidue,k,1),' ',
c     &    gbResidueCoords(countGBResidue,k,2),' ',
c     &    gbResidueCoords(countGBResidue,k,3)
c        call printLine(DEBUG,printstring)
c        k = k + 1
c      enddo


c     save the centre
      do j=1,3
        gbResCentre(countGBResidue,j) = ResCentre(iprot,ires,j)
      enddo
      gbResRadius(countGBResidue) = ResRadius(iprot,ires)
c     Save residue number
      gbResidueList(countGBResidue) = ires
      gbResProtList(countGBResidue) = iprot

      return
      end


