
      integer function findOverride(resnm,ipos)
      include 'dimensions.inc'
      include 'templates.inc'
c########################################################
c
c     This function finds the override template for
c     residue called resnm in chain position ipos
c     or returns 0 if it cannot be found
c
c     (C) Christopher Woods, November 2004
c
c########################################################

      character*(*) resnm
      integer ipos,i
      
      findOverride = 0
      if (ipos.le.0) return
      
      if (NResOverrides(ipos).gt.0) then
        do i=1,NResOverrides(ipos)
          if (ResOverrideNam(ipos,i).eq.resnm) then
            findOverride = ResOverrideIndex(ipos,i)
            return
          endif
        enddo
      endif
      
      return
      end
