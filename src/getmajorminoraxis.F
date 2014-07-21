      subroutine getMajorMinorAxis(ityp,ival1,ival2,itot,imaj,imin,iper)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'templates.inc'
c########################################################
c
c     This routine gets the major/minor/perpendicular
c     axis for a molecule/group (type based on ityp)
c
c     ityp = SOLUTE     isol = ival1
c     ityp = SOLVENT    isvn = ival1
c     ityp = GCSOLUTE   igcs = ival1
c    ***NOT YET IMPLEMENTED*********************************
c    *ityp = PROTEIN    Protein molecule (iprot = ival1)         *
c    *ityp = RESIDUE    Protein residue (iprot=ival1,ires=ival2) *
c    *******************************************************
c
c#########################################################

      integer ityp,ival1,ival2
      integer iprot,ires,isol,isvn,nats,i,j,igcs
      double precision itot(3),imaj(3),imin(3),iper(3)
      double precision ivec(3),posvec(3),negvec(3)
      double precision poslgth2,neglgth2,ideg,itmp,tmplgth,minlgth2,majlgth2
      double precision rad2Deg,veclength,veclength2,dot
      integer nposvec,nnegvec
      
      isol = 0
      isvn = 0
      igcs = 0
      iprot = 0
      ires = 0
      nats = 0

      if (ityp.eq.SOLUTE) then
        isol = ival1
        if (isol.le.0 .or. isol.gt.NSolutes) return
        nats = NSolAtoms(isol)
      else if (ityp.eq.SOLVENT) then
        isvn = ival1
        if (isvn.le.0 .or. isvn.gt.NSolvents) return
        nats = NSvnAtoms(isvn)
      else if (ityp.eq.GCSOLUTE) then
       igcs = ival1
       if (igcs.le.0 .or. igcs.gt.NGCsolutes) return      
       nats = NGCSAtoms(igcs)
      else if (ityp.eq.PROTEIN) then
        iprot = ival1
        if (iprot.le.0 .or. iprot.gt.NProteins) return
        nats = NZAtoms(iprot)
      else if (ityp.eq.RESIDUE) then
        iprot = ival1
        if (iprot.le.0 .or. iprot.gt.NProteins) return
        ires = ival2
        if (ires.le.0 .or. ires.gt.NResidues(iprot)) return
        nats = NSCAtoms(iprot,ires)+3
      endif

      if (nats.le.0) then
        call printLine(WARNING,"Trying to assign dummies to a molecule with no atoms!")
        call printLine(WARNING,"This looks like a problem with your template or coordinate files?")
        return
      else if (nats.lt.2) then
c       assign basic dummies, since there are only 1 (or eventually 2) atoms
        if (nats.eq.1) then
          do 80 i=1,3
            if (ityp.eq.SOLUTE) then
              itot(i) = ZSolCoords(isol,1,i) + 0.3
            else if (ityp.eq.SOLVENT) then
              itot(i) = ZSvnCoords(isvn,1,i)
            else if (ityp.eq.GCSOLUTE) then
              itot(i) = ZGCSCoords(igcs,1,i)
            endif
            imaj(i) = ZERO
            imin(i) = ZERO
            iper(i) = ZERO
80        continue
c         set the major axis up in y, minor along x, and perp along z
          imaj(2) = ONE
          imin(1) = ONE
          iper(3) = ONE
        endif
        return
      endif

c     get the mean geometric centre of the molecule...      
      itot(1) = ZERO
      itot(2) = ZERO
      itot(3) = ZERO
      do 100 i=1,nats
        if (ityp.eq.SOLUTE) then
          do 1011 j=1,3
            itot(j) = itot(j) + ZSolCoords(isol,i,j)
1011      continue
        else if (ityp.eq.SOLVENT) then
          do 1012 j=1,3
            itot(j) = itot(j) + ZSvnCoords(isvn,i,j)
1012      continue
        else if (ityp.eq.GCSOLUTE) then
          do 1013 j=1,3
            itot(j) = itot(j) + ZGCSCoords(igcs,i,j)
1013      continue
        else
          call printLine(WARNING,"Getting the centre of a residue is not yet implemented!")
          call printLine(WARNING,"Maybe you would like to implement it?")
          return
        endif
100   continue

      do 102 i=1,3
        itot(i) = itot(i) / dble(nats)
        posvec(i) = ZERO
        negvec(i) = ZERO
102   continue
      nposvec = 0
      nnegvec = 0
c"      
c     now get the major axis...
c     do this by calculating the vectors of all the atoms from the
c     centre of geometry, and classifying them as 'positive' or
c     'negative' based on whether they have a positive angle with
c     (1,0,0) or a negative angle. (>180<360)
c     The mean average of all of the positive vectors, and the -mean 
c     average of all of the negative vectors should be the mean
c     major axis - the direction can be set as positive towards the 
c     largest mean vector.
      do 200 i=1,nats
c       get the vector
        if (ityp.eq.SOLUTE) then
          do 2011 j=1,3
            ivec(j) = ZSolCoords(isol,i,j) - itot(j)
2011      continue
        else if (ityp.eq.SOLVENT) then
          do 2012 j=1,3
            ivec(j) = ZSvnCoords(isvn,i,j) - itot(j)
2012      continue
         else if (ityp.eq.GCSOLUTE) then
          do 2013 j=1,3
            ivec(j) = ZGCSCoords(igcs,i,j) - itot(j)
2013      continue
        else
          write(printstring,*) "Requesting unknown molecule type! ",ityp 
          call printLine(WARNING,printstring)
          call printLine(WARNING,"This is probably a bug in the program?")
          return
        endif
c"        
c       now get whether it is positive or negative
        tmplgth = veclength(ivec)
        if (tmplgth.gt.ZERO) then
          itmp = dacos(ivec(1)/tmplgth)
          ideg = rad2Deg(itmp)
          if (ideg.gt.THREESIXTY) ideg = ideg - THREESIXTY
          if (ideg.lt.ZERO) ideg = ideg + THREESIXTY
c         positive is -90.0D+00<=ityp<90.0D+00
          if (ideg.ge.-NINETY .and. ideg.lt.NINETY) then
            do 202 j=1,3
              posvec(j) = posvec(j) + ivec(j)
              nposvec = nposvec + 1
202         continue
          else
            do 203 j=1,3
              negvec(j) = negvec(j) + ivec(j)
              nnegvec = nnegvec + 1
203         continue
          endif
        else
          write(printstring,*) "AXIS GENERATION: TMPLGTH <= 0 ",isvn
          call printLine(WARNING,printstring)
        endif
200   continue

      do 204 i=1,3
        if (nposvec.gt.0) posvec(i) = posvec(i) / nposvec
        if (nnegvec.gt.0) negvec(i) = negvec(i) / nnegvec
204   continue

      poslgth2 = veclength2(posvec)
      neglgth2 = veclength2(negvec)
      
      if (poslgth2.gt.neglgth2) then
        if (poslgth2.ne.ZERO) then
          do 210 i=1,3
            imaj(i) = posvec(i)
210       continue
          call vecnormalise(imaj)
        else
           imaj(1) = ONE
           imaj(2) = ZERO
           imaj(3) = ZERO
        endif
      else
        if (neglgth2.ne.ZERO) then
          do 211 i=1,3
            imaj(i) = negvec(i)
211       continue
          call vecnormalise(imaj)
        else
           imaj(1) = ONE
           imaj(2) = ZERO
           imaj(3) = ZERO
        endif
      endif
      
      do 290 i=1,3
        negvec(i) = ZERO
        posvec(i) = ZERO
290   continue
      nnegvec = 0
      nposvec = 0        
c"
c     now do the same, but this time against imaj
      do 300 i=1,nats
        if (ityp.eq.SOLUTE) then
          do 3011 j=1,3
            ivec(j) = ZSolCoords(isol,i,j) - itot(j)
3011      continue
        else if (ityp.eq.SOLVENT) then
          do 3012 j=1,3
            ivec(j) = ZSvnCoords(isvn,i,j) - itot(j)
3012      continue
        else if (ityp.eq.GCSOLUTE) then
          do 3013 j=1,3
            ivec(j) = ZGCSCoords(igcs,i,j) - itot(j)
3013      continue
        else
          call printLine(WARNING,"Getting the centre of a residue is not yet implemented!")
          call printLine(WARNING,"Maybe you would like to implement it?")
          return
        endif
c" 
        tmplgth = veclength(ivec)
        if (tmplgth.gt.ZERO) then
          ideg = tmplgth*veclength(imaj)
          itmp = dot(imaj,ivec)
          itmp = dacos(itmp)
          itmp = rad2Deg(itmp)
          ideg = itmp / ideg
          if (ideg.gt.THREESIXTY) ideg = ideg - THREESIXTY
          if (ideg.lt.ZERO) ideg = ideg + THREESIXTY
c         positive is -90.0D+00<=ityp<90.0D+00
          if (ideg.ge.-NINETY .and. ideg.lt.NINETY) then
            do 302 j=1,3
              posvec(j) = posvec(j) + ivec(j)
              nposvec = nposvec + 1
302         continue
          else
            do 303 j=1,3
              negvec(j) = negvec(j) + ivec(j)
              nnegvec = nnegvec + 1
303         continue
          endif
        else
          write(printstring,*) "AXIS GENERATION: ZERO TMPLGTH ",isvn
          call printLine(INFO,printstring)
          call printLine(WARNING,printstring)
        endif
300   continue

      do 304 i=1,3
        if (nposvec.gt.0) posvec(i) = posvec(i) / nposvec
        if (nnegvec.gt.0) negvec(i) = negvec(i) / nnegvec
304   continue

      poslgth2 = veclength2(posvec)
      neglgth2 = veclength2(negvec)

      if (poslgth2.gt.neglgth2) then
        if (poslgth2.gt.TINY) then
          do 310 i=1,3
            imin(i) = posvec(i)
310       continue
          call vecnormalise(imin)
        else
           imin(1) = ZERO
           imin(2) = ONE
           imin(3) = ZERO
        endif
      else
        if (neglgth2.gt.TINY) then
          do 311 i=1,3
            imin(i) = negvec(i)
311       continue
          call vecnormalise(imin)
        else
           imin(1) = ZERO
           imin(2) = ONE
           imin(3) = ZERO
        endif
      endif

c     now get the lengths of the two vectors...
      majlgth2 = veclength2(imaj)
      minlgth2 = veclength2(imin)
      if (majlgth2.lt.TINY) then
c       zero vector!
c       set equal to 1,0,0
        imaj(1) = ONE
        imaj(2) = ZERO
        imaj(3) = ZERO
      else if (minlgth2.lt.TINY) then
c       zero vector!
c       set equal to 0,1,0
        imin(1) = ZERO
        imin(2) = ONE
        imin(3) = ZERO
      endif
      
      if (minlgth2.gt.majlgth2) then
c       we need to swap them over!
        do 340 i=1,3
          ivec(i) = imaj(i)
          imaj(i) = imin(i)
          imin(i) = ivec(i)
340     continue
      endif

      do 390 i=1,3
        ivec(i) = imin(i)
390   continue

      call vecnormalise(ivec)
      call vecnormalise(imaj)
      
      call getAxisVectors(imaj,ivec,imin,iper)
      
c     for now, check that the axis are all perpendicular
      if (abs(dot(imaj,imin)).gt.SMALL .or. abs(dot(imaj,iper)).gt.SMALL
     .                                 .or. abs(dot(imin,iper)).gt.SMALL) then
        call printLine(WARNING,"Something went very wrong with axis generation!")
        call printLine(WARNING,"This is a program bug, or numerical instability problem!")
        return
      endif

      return
      end
