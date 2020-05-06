      subroutine markgbupdate
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'gb.inc'
c######################################
c     Determine which part of the system
c     needs to have its gb energy updated
c     
c     Julien Michel November 2004
c######################################
      integer i,j,k,nres,resstart,resend,nats

      if (NProteins .gt. 0) then
        do 10 i=1,NProteins
          nres = NResidues(i)
          do 12 j=1,nres
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)
c           Assume no need to redo it's GB energy
            calcResidueGB(i,j) = .false.
            do 14 k=resstart,resend
              ProtAtomdoBorn(i,k) = .false.
              if ( MovedResidue(i,j) .gt. 0
     &          .or. 
     &          (abs(ZProtEffBorn(i,k)-OldZProtEffBorn(i,k)) .gt. gbthreshold ) 
     &          .or.
     &          (abs(ZProtEffBornF(i,k)-OldZProtEffBornF(i,k)) .gt. gbthreshold )
     &          .or.
     &          (abs(ZProtEffBornB(i,k)-OldZProtEffBornB(i,k)) .gt. gbthreshold) ) 
     &           then
c               write(printstring,*)'Need update for atom ',k,
c               &            ' in residue',j,' protein ',i
c                call printLine(DEBUG,printstring)
                ProtAtomdoBorn(i,k) = .true.
                calcResidueGB(i,j) = .true.
              endif              
 14         continue
 12       continue
 10     continue
      endif

      if (NSolutes.gt. 0) then
        do 20 i=1,NSolutes
          nats = NSolAtoms(i)
c         Assume no need to redo GB
          calcSoluteGB(i) = .false.
          do 30 j=1,nats
            SolAtomdoBorn(i,j) = .false.
c            write(printstring,*) 'Zsol ',ZSolEffBorn(i,j),
c     .        ' Old ',OldZSolEffBorn(i,j),
c     .        ' diff ',abs( ZSolEffBorn(i,j)-OldZSolEffBorn(i,j))
c            call printLine(DEBUG,printstring)
            if ( MovedSolute(i) .gt. 0 .or. 
     &        abs( ZSolEffBorn(i,j)-OldZSolEffBorn(i,j)) .gt. gbthreshold .or.
     &        abs( ZSolEffBornF(i,j)-OldZSolEffBornF(i,j)) .gt. gbthreshold .or.
     &        abs( ZSolEffBornB(i,j)-OldZSolEffBornB(i,j)) .gt. gbthreshold ) then
              SolAtomdoBorn(i,j) = .true.
              calcSoluteGB(i) = .true.
            endif
 30       continue
 20     continue
      endif
      
      if (NSolvents .gt. 0) then
        do 40 i=1,NSolvents
          nats = NSvnAtoms(i)
c         Assume no need to redo it's GB energy
          calcSolventGB(i) = .false.
          do 50 j=1,nats
            SvnAtomdoBorn(i,j) = .false.
            if (MovedSolvent(i) .gt. 0 .or.
     &        abs( ZSvnEffBorn(i,j)-OldZSvnEffBorn(i,j)) .gt. gbthreshold .or.
     &        abs( ZSvnEffBornF(i,j)-OldZSvnEffBornF(i,j)) .gt. gbthreshold .or.
     &        abs( ZSvnEffBornB(i,j)-OldZSvnEffBornB(i,j)) .gt. gbthreshold ) then
              SvnAtomdoBorn(i,j) = .true.
              calcSolventGB(i) = .true.
            endif
 50       continue
 40     continue
      endif

      
c      do i=1,NProteins
c        do j=1,NResidues(i)
c          do k=BBZid(i,j,1),BBZid(i,j,4)
c            write(printstring,*) 'Prot ',i,' Res ',j,' at ',k,
c     .        ' update ',ProtAtomdoBorn(i,k)
c            call printLine(DEBUG,printstring) 
c          enddo
c        enddo
c      enddo
c      do i=1,NSolutes
c        do j=1,NSolAtoms(i)
c          write(printstring,*) 'Sol ',i,' at ',j,' update ',SolAtomdoBorn(i,j)
c          call printLine(DEBUG,printstring)
c        enddo
c      enddo
c
c      do i=1,NSolvents
c        do j=1,NSvnAtoms(i)
c          write(printstring,*) 'Svn ',i,' at ',j,' update ',SvnAtomdoBorn(i,j)
c          call printLine(DEBUG,printstring)
c        enddo
c      enddo


      return
      end
