      integer function findTemplateDihedral(at1,at2,at3,at4,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#####################################################
c
c     This routine searchs for a dihedral angle between amber
c     atom names at1-at2-at3-at4 and returns the index 
c     number in the array - it returns 0 if the 
c     angle does not exist...
c
c     This function assumes that the angle array
c     has atoms internally sorted...
c
c     (C) Christopher Woods  23 January 2002
c
c######################################################

      character*4 at1,at2,at3,at4
      integer i,itmpl

      findTemplateDihedral = 0
      if (NTemplateDihedrals(itmpl).le.0) return
      
      call sortTemplateDihedral(at1,at2,at3,at4)
      
      do i=1,NTemplateDihedrals(itmpl)
        if (at1.eq.TemplateDihName(itmpl,i,1)) then
          if (at2.eq.TemplateDihName(itmpl,i,2)) then
            if (at3.eq.TemplateDihName(itmpl,i,3)) then
              if (at4.eq.TemplateDihName(itmpl,i,4)) then
                findTemplateDihedral = i
                return
              endif
            endif
          endif
        endif
      enddo

      return
      end
