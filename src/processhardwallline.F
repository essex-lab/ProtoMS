      subroutine processHardwallLine(check)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'labels.inc'
      include 'extraenergies.inc'
c#########################################
c
c     Process the 'hardwall' input line - if check
c     is true then only print out what we would do.
c
c     The line has been split into words that are held
c     in the 'words' array
c
c     The ID must have been set previously for this to work
c
c
c     hardwall idnum ixcen iycen izcen irad ibigeng
c
c     --> Add a hardwall potential to idnum=#, this is set by the add command and defines 
c        type, atom number, residue number etc.
c        so before using the 'hardwall' line need: 'id add 1 SOLUTE 1 C19 3T4' (for your atom of interest)
c         Potential is centered at ixcen iycen izcen with radius irad 
c         Atoms within irad experience ibigeng
c         Both this subroutine and addid need to be called in the chunk section
c
c        For now this only works for solutes
c
c     IMPORTANT
c     Once a restraint is added, it cannot be cleared
c     This is because it is difficult to remove them cleanly 
c     and get the correct running averages and total energies
c     Also, you may not get what you want if you start to simulate before
c     adding the restraint ...
c
c      J. Taylor July 2005
c       Adapted from: processrestraintline.F and addcartesianharmonicrestraint.F 
c       By: Julien Michel February 2005
c     Modification: 
c       By A.I.Cabedo Martinez 2014
c       The input can now be hardwall idnum irad ibigeng
c       When this is the case, the current position of the geometry centre of the molecule of interest
c       is assigned to ixcen iycen izcen
c
c
c#########################################
      logical check,allowprotein
      integer idnum
      integer getInt
      double precision getDouble
      double precision ixcen,iycen,izcen,irad,ibigeng
      integer zid1
      integer typeid1
      integer molid1
      character*4 resnameid1
      integer resnumid1
      character*4 nameid1
      
      if (nwords.ne.7.and.nwords.ne.8.and.nwords.ne.4.and.nwords.ne.5) then
        call printLine(WARNING,"Hardwall line not correctly implemented.")
        call join(1,nwords,printstring)
      call printLine(WARNING,"Correct useage: hardwall idnum ixcen iycen izcen irad ibigeng (allowprotein)",
     .                       " or hardwall idnum irad ibigeng (allowprotein)")
        call printLine(WARNING,printstring)
        return
      endif
c         
c Get all the paramters that are always present  
      if (nwords.ge.7) then
        idnum = getInt(words(2),0)
        ixcen = getDouble(words(3),0)
        iycen = getDouble(words(4),0)
        izcen = getDouble(words(5),0)
        irad = getDouble(words(6),0)
        ibigeng = getDouble(words(7),0)
      elseif (nwords.ge.4) then

        idnum = getInt(words(2),0)
        irad = getDouble(words(3),0)
        ibigeng = getDouble(words(4),0)
      endif

c Set allowprotein to true if option is present
      allowprotein = .false.
      if (nwords.eq.8) then
         if (words(8).eq.'allowprotein') then
            allowprotein = .true.
         endif 
      else if (nwords.eq.5) then
         if (words(5).eq.'allowprotein') then
            allowprotein = .true.
         endif 
      endif


c Set cubic hardwall to true if option is present
      cubichardwall=.false.
      if (nwords.eq.5) then
         if (words(5).eq.'cubichardwall') then
            cubichardwall = .true.
         endif 
      endif

c       Make sure the IDs are valid
        if ( idnum .ne. 0 .and. LabelID(idnum) .eq. 0) then
          write(printstring,*)' No ID = ',idnum,' has been set ! ',
     .      ' cannot process the following line'
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)                 
          return
        else if (idnum .eq. 0) then
          continue
        else
          zid1 = LabelID(idnum)
          typeid1 = LabelType(idnum) 
          molid1 = LabelMolecule(idnum) 
          resnameid1 = LabelResidueName(idnum) 
          resnumid1 = LabelResidueNumber(idnum)
          nameid1 = LabelName(idnum)
        endif

      if (nwords.eq.4.or.nwords.eq.5) then
        ixcen = SolGeomCentre(molid1,1)
        iycen = SolGeomCentre(molid1,2)
        izcen = SolGeomCentre(molid1,3)
      endif

      if (typeid1 .eq. PROTEIN) then
        call printLine(WARNING,"Hardwall potentials can not be applied to proteins at this time")
        return
      elseif (typeid1 .eq. SOLVENT) then
        call printLine(WARNING,"Hardwall potentials can not be applied to solvents at this time")
        return
      elseif (typeid1 .eq. SOLUTE) then
            write(printstring,*)'Adding hardwall potential of ',ibigeng,' kcal/mol on atom with id ',idnum,' with center ',
     . ixcen,iycen,izcen,' and radius ',irad
            call printLine(INFO,printstring)
      endif
c          

      nHardwallPotential = nHardwallPotential + 1
      HardwallCenter(nHardwallPotential,1) = ixcen
      HardwallCenter(nHardwallPotential,2) = iycen
      HardwallCenter(nHardwallPotential,3) = izcen
      HardwallRadius(nHardwallPotential) = irad
      HardwallMolecule(nHardwallPotential) = molid1
      HardwallextraEnergy(nHardwallPotential) = ibigeng
      HardwallAllowProtein(nHardwallPotential) = allowprotein
      
           return
            end
      
