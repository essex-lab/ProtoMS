      subroutine processRestraintLine(check)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'labels.inc'
c#########################################
c
c     Process the 'restrain' input line - if check
c     is true then only print out what we would do.
c
c     The line has been split into words that are held
c     in the 'words' array
c
c     The ID must have been set previously for this to work
c
c
c     restraint add 1 cartesian harmonic X,Y,Z K
c     --> Add a parabolic restraint that acts on the cartesian coordinates 
c     of atom ID 1 and is centered on cartesian coordinates X,Y,Z and with 
c     force constant k 
c
c     restraint add 1 dihedral harmonic Rad,K
c     --> add a parabolic restraint centered on Rad with force constant K
c     that acts on the zmatrix torsion of atom with ID 1
c     For this to work, the atom must have a zmatrix torsion
c
c     restraint add 1-2-3-4 dihedral harmonic Rad,K
c     --> Add a parabolic restraint centered on Rad, with force constant K
c     that acts on the torsion defined by the four atoms wich have ID 1,2,3,4
c     This is a more general way of defining restraints which are not part 
c     of a zmatrix
c
c     The same thing apply to 
c     1 bond 
c     1-5 bond 
c     7 angle 
c     4-3-1 angle 
c     etc...
c
c     IMPORTANT
c     Once a restraint is added, it cannot be cleared
c     This is because it is difficult to remove them cleanly 
c     and get the correct running averages and total energies
c     Also, you may not get what you want if you start to simulate before
c     adding the restraint ...
c
c     Julien Michel February 2005
c
c#########################################
      logical check
      integer id1,id2,id3,id4
      integer getInt
      double precision getDouble
      double precision xrestr,yrestr,zrestr,krestr,bndrestr,dihrestr,wrestr
      integer zid1,zid2,zid3,zid4
      integer typeid1,typeid2,typeid3,typeid4
      integer molid1,molid2,molid3,molid4
      character*4 resnameid1,resnameid2,resnameid3,resnameid4
      integer resnumid1,resnumid2,resnumid3,resnumid4
      character*4 nameid1,nameid2,nameid3,nameid4
      integer nids
      double precision deg2rad
      logical getLogical, scale

      if (check) then
         write(printstring,*) "Processing restraint lines in checking mode is not supported."
         call printLine(WARNING,printstring)
         return
      endif
      
      if (nwords.lt.4) then
        call printLine(WARNING,"Could not interpret restraint line")
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return
      endif      
      call lowercase(words(2))
      if (words(2) .eq. 'add') then
c       We use split to get the number of atoms involved in this restraint
c       note that if you pass more atoms than what the restrain can handle,
c       then these atoms will be ignored
        call split2(words(3),"-")
        id1 = getInt(words2(1),0)
        id2 = getInt(words2(2),0)
        id3 = getInt(words2(3),0)
        id4 = getInt(words2(4),0)
        nids = nwords2
c        write(printstring,*) 'id1 ',id1,' id2 ',id2,' id3 ',id3,' id4 ',id4
c        call printLine(DEBUG,printstring)
c       Make sure the IDs are valid
        if ( id1 .ne. 0 .and. LabelID(id1) .eq. 0) then
          write(printstring,*)' No ID = ',id1,' has been set ! ',
     .      ' cannot process the following line'
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)                 
          return
        else if (id1 .eq. 0) then
          continue
        else
          zid1 = LabelID(id1)
          typeid1= LabelType(id1) 
          molid1 = LabelMolecule(id1) 
          resnameid1 = LabelResidueName(id1) 
          resnumid1 = LabelResidueNumber(id1)
          nameid1 = LabelName(id1)
        endif

        if ( id2 .ne. 0 .and. LabelID(id2) .eq. 0) then
          write(printstring,*)' No ID = ',id2,' has been set ! ',
     .      ' cannot process the following line'
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)                 
          return
        else if (id2 .eq. 0) then
          continue
        else
          zid2 = LabelID(id2)
          typeid2= LabelType(id2) 
          molid2 = LabelMolecule(id2) 
          resnameid2 = LabelResidueName(id2) 
          resnumid2 = LabelResidueNumber(id2)
          nameid2 = LabelName(id2)
        endif

        if ( id3 .ne. 0 .and. LabelID(id3) .eq. 0) then
          write(printstring,*)' No ID = ',id3,' has been set ! ',
     .      ' cannot process the following line'
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)                 
          return
        else if (id3 .eq. 0) then
          continue
        else
          zid3 = LabelID(id3)
          typeid3= LabelType(id3) 
          molid3 = LabelMolecule(id3) 
          resnameid3 = LabelResidueName(id3)
          resnumid3 = LabelResidueNumber(id3)
          nameid3 = LabelName(id3)
        endif 
       
        if ( id4 .ne. 0 .and. LabelID(id4) .eq. 0) then
          write(printstring,*)' No ID = ',id4,' has been set ! ',
     .      ' cannot process the following line'
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)                 
          return
        else if (id4 .eq. 0) then
          continue
        else
          zid4 = LabelID(id4)
          typeid4= LabelType(id4) 
          molid4 = LabelMolecule(id4) 
          resnameid4 = LabelResidueName(id4)
          resnumid4 = LabelResidueNumber(id4)
          nameid4 = LabelName(id4)
        endif     
     
        call lowercase(words(4))
        if (words(4) .eq. 'cartesian') then 
c
c         Cartesian restraints apply on only one atom , so any id2,id3,id4
c         will simply be ignored
c
          call lowercase(words(5))          
          if ( words(5) .eq. 'harmonic') then
            xrestr = getDouble(words(6),0)
            yrestr = getDouble(words(7),0)
            zrestr = getDouble(words(8),0)
            krestr = getDouble(words(9),0)
            if ( words(10) .eq. 'lambda' ) then
               scale = .true.
            else
               scale = .false.
            endif
            write(printstring,10) words(4),words(5),id1,xrestr,yrestr,zrestr,krestr,scale
            call printLine(INFO,printstring)
 10         Format("Adding ",A12," ",A8," restraint on ID ",I4,
     .        " center ",3(F6.3,1X)," K ",F6.2," lambda scaling =",L)
            call addCartesianHarmonicRestraint(typeid1,molid1,zid1,
     .        xrestr,yrestr,zrestr,krestr,scale)
          else if (words(5) .eq. 'flatbottom') then
            xrestr = getDouble(words(6),0)
            yrestr = getDouble(words(7),0)
            zrestr = getDouble(words(8),0)
            krestr = getDouble(words(9),0)
            wrestr = getDouble(words(10),0)
            write(printstring,11) words(4),words(5),id1,xrestr,yrestr,zrestr,krestr,wrestr
            call printLine(INFO,printstring)
 11         Format("Adding ",A12," ",A8," restraint on ID ",I4,
     .        " center ",3(F6.3,1X)," K ",F6.2," width ",F6.2)
            call addCartesianFlatBottomRestraint(typeid1,molid1,zid1,
     .        xrestr,yrestr,zrestr,krestr,wrestr)
          else 
            write(printstring,*) 'Could not interpret functional ',
     .        ' form of restraint'
            call join(1,nwords,printstring)
            call printLine(WARNING,printstring)                 
            return
          endif
        else if (words(4) .eq. 'bond') then
          call lowercase(words(5))          
          if ( words(5) .eq. 'harmonic') then
            krestr = getDouble(words(6),0)
            bndrestr = getDouble(words(7),0)
            write(printstring,21) words(4),words(5),id1,id2,krestr,bndrestr
            call printLine(INFO,printstring)
 21         Format("Adding ",A12," ",A8," restraint on ID ",I4,"-",I4,
     .        " force ",F6.3," equilibrium ",F6.3)
            call addBondHarmonicRestraint(typeid1,molid1,zid1,
     .        typeid2,molid2,zid2,krestr,bndrestr)
          else 
            write(printstring,*) 'Could not interpret functional ',
     .        ' form of restraint'
            call join(1,nwords,printstring)
            call printLine(WARNING,printstring)                 
            return            
          endif
        else if (words(4) .eq. 'angle') then
          write(printstring,*) 'Sorry, angle restraints not in yet'
          call printLine(WARNING,printstring)  
        else if (words(4) .eq. 'dihedral') then
          call lowercase(words(5))
          if ( words(5) .eq. 'harmonic') then
            dihrestr = getDouble(words(6),0)
            dihrestr = deg2rad(dihrestr)
            krestr = getDouble(words(7),0)
            write(printstring,12) words(4),words(5),id1,id2,id3,id4,dihrestr,krestr
            call printLine(INFO,printstring)
 12         Format("Adding ",A12," ",A8," restraint on ID ",4(I4,1X),
     .        " equil ",F6.3," rad  K ",F6.2)
            if (id2 .eq.0 .or. id3 .eq. 0 .or. id4 .eq. 0) then
              if (typeid1 .eq. SOLVENT) then
                write(printstring,*) 'Cannot set a dihedral restraint '
     .            ,' for a SOLVENT from a zmatrix ! '
                call printLine(WARNING,printstring)
                return
              else if (typeid1 .eq. PROTEIN) then
                if (nameid1 .eq. 'n' 
     .            .or. nameid1 .eq. 'c' 
     .            .or. nameid1 .eq. 'ca'
     .            .or. nameid1 .eq. 'o') then
                  write(printstring,*) 'Cannot set a dihedral restraint '
     .              ,' for a backbone PROTEIN atom from its zmatrix ! '
                  call printLine(WARNING,printstring)
                  return                
                endif
              endif
              call addDihedralHarmonicRestraint(typeid1,molid1,resnumid1,zid1,
     .          0,0,0,0,0,0,0,0,0,0,0,0,
     .          dihrestr,krestr)
            else 
              call addDihedralHarmonicRestraint(typeid1,molid1,resnumid1,zid1,
     .          typeid2,molid2,resnumid2,zid2,typeid3,molid3,resnumid3,zid3,
     .          typeid4,molid4,resnumid4,zid4,
     .          dihrestr,krestr)              
            endif
           
          endif
        else   
          write(printstring,*) 'Could not interpret type of restraint'
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)          
          return
        endif
      else 
        write(printstring,*) 'Could not interpret restraint line'
        call printLine(WARNING,printstring)
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)          
        return        
      endif

      return
      end
      
