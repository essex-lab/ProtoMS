      character*2 function elementName(iel)
      implicit none
c##############################################################
c
c     This function will return the element name of the element
c     with proton number 'iel'
c
c     (C) Christopher Woods 20 June 2002
c
c##############################################################

      integer iel
      character*2 elnames(20)
      
c     data arranged by periodic table z number - only goes up to Ca,
c     anything outside range is a dummy (DM)
      data elnames / 'H ',                                   'He',
     .               'Li','Be',     'B ','C ','N ','O ','F ','Ne',
     .               'Na','Mg',     'Al','Si','P ','S ','Cl','Ar',
     .               'K ','Ca' /
       
      if (iel.le.0 .or. iel.gt.20) then
        elementName = 'DM'
        return
      endif
       
      elementName = elnames(iel)
      return
      end
