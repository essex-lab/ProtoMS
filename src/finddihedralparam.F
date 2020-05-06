      integer function findDihedralParam(at1,at2,at3,at4)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#########################################################
c
c     This function returns the dihedral parameter
c     for amber atoms at1-at2-at3-at4, or 0 if
c     this is a dummy dihedral or it is not in
c     the dihedral parameter database
c
c     (C) Christopher Woods
c
c#########################################################

      character*2 at1,at2,at3,at4,tat1,tat2,tat3,tat4
      integer findDihedralIndex,idih
      character *2 x1,x2
      findDihedralParam = findDihedralIndex(at1,at2,at3,at4)
c        write(printstring,*) 'TORSION for ',at1,at2,at3,at4,
c     .    ' = ',findDihedralParam
c        call printLine(WARNING,printstring)

      if (findDihedralParam .le. 0 .and. at1.ne.'dm' .and. at2.ne.'dm' .and.
     .  at3.ne.'dm' .and. at4.ne.'dm' .and. at1.ne.'DU' .and. 
     .  at2.ne.'DU' .and. at3.ne.'DU' .and. at4.ne.'DU') then
c       Not found...Can we find a matching generic torsion ? 
c       Do this only if the torsion doesn't involve dummy atoms 
        x1 = 'X '
        tat1 = at2
        tat2 = at3
        x2 = 'X '
        findDihedralParam = findDihedralIndex(x1,tat1,tat2,x2)
c       write(printstring,*) 'GENERIC TORSION for ',at1,at2,at3,at4,
c     .    ' = ',findDihedralParam
c       call printLine(WARNING,printstring)
      endif
      
c     The dihedral was not found! - this is ok if either of the dihedral
c     atoms are dummies
      if (findDihedralParam.le.0 .and. at1.ne.'dm' .and. at2.ne.'dm' .and.
     .  at3.ne.'dm' .and. at4.ne.'dm' .and. at1.ne.'DU' .and. 
     .  at2.ne.'DU' .and. at3.ne.'DU' .and. at4.ne.'DU') then
        write(printstring,10) at1,at2,at3,at4
10      format("Could not find dihedral parameter between amber atoms ",a2,"-",a2,"-",a2,"-",a2,
     .           ". Will treat this dihedral as a dummy. This dihedral should not be sampled!")
        call printLine(WARNING,printstring)
c       now add this dihedral as a dummy to prevent future warnings
        idih = NDihParam + 1
        
        if (idih.le.MAXDIHPARAM) then
          tat1 = at1
          tat2 = at2 
          tat3 = at3
          tat4 = at4
          call sortDihedral(tat1,tat2,tat3,tat4)
          
          DihName(idih,1) = tat1
          DihName(idih,2) = tat2
          DihName(idih,3) = tat3
          DihName(idih,4) = tat4
          DihParamIndex(idih) = 0
          NDihParam = idih
        endif
      endif
            
      if (findDihedralParam.gt.0) 
     .    findDihedralParam = DihParamIndex(findDihedralParam)
     
      return
      end
