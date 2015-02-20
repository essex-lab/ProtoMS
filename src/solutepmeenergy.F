      subroutine solutePmeEnergy(isol,exclnrg,selfnrg,charges,x,y,z,natoms,maxpmeatoms)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'flags.inc'
      include 'quickpot.inc'
c#######################################################################
c
c     This subroutine calculates necessary energies for 
c     particle-mesh Ewald for solute isol
c
c     (C) Samuel Genheden, October 2013
c
c#######################################################################

      integer isol,natoms,maxpmeatoms
      double precision exclnrg,q2sum,selfnrg,exclsolnrg,selfsolnrg
      double precision charges(maxpmeatoms),x(maxpmeatoms),y(maxpmeatoms),z(maxpmeatoms)
      integer i,j,k,innb
      include 'nbenergyvar2.inc'
      double precision charge,mindist      
      double precision scl,sclf,sclb
      logical take

      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
      mindist = 99999.0
      softcore = .false.
      ddd = .false.

      exclsolnrg = ZERO
      selfsolnrg = ZERO
      q2sum = ZERO

      if (IgnoreSolute(isol)) return
      
c     number non-bonded pairs in the solute molecule
      innb = NSolNB(isol)

c     no periodic boundaries in an intramolecular NB calculation
      boundx = .false.
      boundy = .false.
      boundz = .false.

c     Calculate exclusion energy and sum of charges squared
      funcisol1 = isol
      funcisol2 = isol
      do i=1,NSolAtoms(isol)
        funcat1 = i
        do j=i+1,NSolAtoms(isol)
          funcat2 = j
          
c         Check to see if this pair is among the non-bonded pairs          
          take = .true.
          do k=1,innb
c            if (i.eq.8) write(6,*)"NB",SolNB(isol,k,1),SolNB(isol,k,2),i,j
            if (((SolNB(isol,k,1).eq.i).and.(SolNB(isol,k,2).eq.j)).or.((SolNB(isol,k,2).eq.i).and.(SolNB(isol,k,1).eq.j)))
     .         take = .false.
            if (SolNB14(isol,k)) take = .true.
          enddo

          if (take) then
c            write(6,*)'E',i,j
c           initialise the calculation
            include 'nbinitial.inc'
c           set up for a solute-solute interaction
            include 'nbsolsol.inc'
c           set up the charge
            include 'qinitial.inc'
c           actually perform the calculation
            include 'qenergy.inc'
            exclsolnrg = exclsolnrg + funccnrg
          endif
        enddo


c       Now calculate the sum of the charges squared (for self-energy)
        ipar1 = ZSolPar(funcisol1,funcat1)
        ipar2 = ZSolPar(funcisol2,funcat2)
c       set up the charge
        include 'qinitial.inc'
        q2sum = q2sum + funcichg1*funcichg1
       
c       Load arrays for mesh calculation
        natoms = natoms + 1
        x(natoms) = ZSolCoords(isol,i,1)/10.0d0
        y(natoms) = ZSolCoords(isol,i,2)/10.0d0
        z(natoms) = ZSolCoords(isol,i,3)/10.0d0
        charges(natoms) = funcichg1
      enddo

      selfsolnrg = -q2sum*FAC4PIEPS0*pme_alpha/sqrt(PI)
      exclsolnrg = -exclsolnrg

      selfnrg = selfnrg + selfsolnrg
      exclnrg = exclnrg + exclsolnrg
      return
      end
