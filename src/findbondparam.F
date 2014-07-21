      integer function findBondParam(at1,at2)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'enums.inc'
c######################################################
c
c     This function returns the parameter number
c     for the bond between amber atoms at1 and at2,
c     or 0 if this is a dummy bond or it doesn't
c     exist in the database
c
c     (C) Christopher Woods, October 2004
c
c######################################################

      character*2 at1,at2,tat1,tat2
      integer findBondIndex,ibnd
      
      findBondParam = findBondIndex(at1,at2)
      

c     The bond was not found! - this is ok if either of the bond
c     atoms are dummies
      if (findBondParam.le.0 .and. at1.ne.'dm' .and. at2.ne.'dm') then
        write(printstring,10) at1,at2
10      format("Could not find bond parameter between amber atoms ",a2,"-",a2,
     .           ". Will treat this bond as a dummy. This bond should not be sampled!")
        call printLine(WARNING,printstring)
        
c       now add the bond as a dummy so that we don't get flooded with warnings!
        ibnd = NBndParam + 1
        if (ibnd.le.MAXBNDPARAM) then
          tat1 = at1
          tat2 = at2
          call sortBond(tat1,tat2)
          
          BndName(ibnd,1) = tat1
          BndName(ibnd,2) = tat2
          BndParamIndex(ibnd) = 0
          NBndParam = ibnd
        endif
      endif
      
      if (findBondParam.gt.0) findBondParam = BndParamIndex(findBondParam)
      
      return
      end
