      integer function findAngleParam(at1,at2,at3)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'enums.inc'
c######################################################
c
c     This function returns the angle parameter
c     number for the angle between amber atoms
c     at1,at2,at3, or zero if this is a dummy
c     angle or this angle does not exist
c
c     (C) Christopher Woods, October 2004
c
c######################################################

      character*2 at1,at2,at3
      character*2 tat1,tat2,tat3
      integer findAngleIndex,iang
      
      findAngleParam = findAngleIndex(at1,at2,at3)
      
c     The angle was not found! - this is ok if either of the angle
c     atoms are dummies
      if (findAngleParam.le.0 .and. at1.ne.'dm' .and. at2.ne.'dm' .and. at3.ne.'dm') then
        write(printstring,10) at1,at2,at3
10      format("Could not find angle parameter between amber atoms ",a2,"-",a2,"-",a2,
     .           ". Will treat this angle as a dummy. This angle should not be sampled!")
        call printLine(WARNING,printstring)
        
c       now add this angle as a dummy to prevent future warnings
        iang = NAngParam + 1
        
        if (iang.le.MAXANGPARAM) then
          tat1 = at1
          tat2 = at2 
          tat3 = at3
          call sortAngle(tat1,tat2,tat3)
         
          AngName(iang,1) = tat1
          AngName(iang,2) = tat2
          AngName(iang,3) = tat3
          AngParamIndex(iang) = 0
          NAngParam = iang
        endif
      endif
 
      if (findAngleParam.gt.0) findAngleParam = AngParamIndex(findAngleParam)
      
      return
      end
