      subroutine getSoluteNonBonded(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'enums.inc'
c#########################################################
c
c     This routine gets all of the solute non-bonded 
c     pairs, and highlights all of the 1-4 
c     interactions.
c
c     (C) Christopher Woods
c
c#########################################################

      integer isol,i,j,k,l,nats,nbnds
      logical bonded(MAXSOLUTEATOMS,MAXSOLUTEATOMS)
      integer bondtyp,iang,idih
      
      NSolNB(isol) = 0

      nats = NSolAtoms(isol)
      if (nats.le.2) return

      write(printstring,"(a,i4)") "Assigning non-bonded atoms for solute ",isol
      call printLine(INFO,printstring)
     
c     initialise the bonded array
      do i=1,nats
        do j=1,nats
          bonded(i,j) = .false.
        enddo
        bonded(i,i) = .true.
      enddo
      
c     only atoms that are part of explicit, non-dummy bonds are bonded together
      nbnds = NSolBonds(isol)
      if (nbnds.gt.0) then
        do i=1,nbnds
          if (.not.SolBondDummy(isol,i)) then
            bonded(SolBond(isol,i,1),SolBond(isol,i,2)) = .true.
            bonded(SolBond(isol,i,2),SolBond(isol,i,1)) = .true.
          endif
        enddo
      endif
           
c     now loop through the solute atoms and work out what is bonded/angled/dihedraled to what...       
      do i=1,nats-1
        do j=i+1,nats
          bondtyp = 0
          
          if (bonded(i,j)) then
c           the atoms are bonded            
            bondtyp = 1
            goto 100
          else
c           these two atoms are not bonded together! - see if they are angled
            do k=1,nats
              if (j.ne.k .and. i.ne.k) then
                if (bonded(i,k).and.bonded(k,j)) then
c                 these two atoms are angled! (i-k-j) - add these as an implicit angle
                  bondtyp = 2
                  iang = k
                  goto 100
                endif
              endif
            enddo

            do k=1,nats
              if (j.ne.k .and. i.ne.k) then
c               see if these two atoms are dihedraled
                if (bonded(i,k)) then
                  do l=1,nats
                    if (i.ne.l .and. j.ne.l .and. k.ne.l) then
                      if (bonded(k,l).and.bonded(l,j)) then                
c                       these atoms form a dihedral (i-k-l-j)
                        bondtyp = 3
                        iang = k
                        idih = l
                        goto 100
                      endif
                    endif
                  enddo
                endif
              endif
            enddo
            
          endif
c         break to here if find that the atoms are in an bond/angle/dihedral          
100       continue        
          if (bondtyp.eq.1) then
c           the atoms are bonded
          else if (bondtyp.eq.2) then          
c           the atoms are angled
            call addSoluteAngle(isol,i,iang,j,.false.,.false.,ZERO,-1,-1)
          else if (bondtyp.eq.3) then
c           the atoms are dihedraled
            call addSoluteDihedral(isol,i,iang,idih,j,.false.,.false.,ZERO,-1,-1)
            call addSoluteNB(isol,i,j,.true.)          
          else        
c           the atoms are not bonded
            call addSoluteNB(isol,i,j,.false.)
          endif
        
        enddo
      enddo

      return
      end
