      integer function findUreyBradleyParam(at1,at2,at3)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'enums.inc'
c#########################################################
c
c     This returns the uby parameter number of the uby
c     between amber atoms at1-at2-at3, or 0 if this 
c     is a dummy uby, or it isn't in the uby database
c
c     (C) Christopher Woods, October 2004
c
c#########################################################

      character*2 at1,at2,at3,tat1,tat2,tat3
      integer findUreyBradleyIndex,iuby
      
      findUreyBradleyParam = findUreyBradleyIndex(at1,at2,at3)
      
c     The uby was not found! - this is ok if either of the uby
c     atoms are dummies
      if (findUreyBradleyParam.le.0 .and. at1.ne.'dm' .and. at2.ne.'dm' .and. at3.ne.'dm') then
        write(printstring,10) at1,at2,at3
10      format("Could not find Urey-Bradley parameter between amber atoms ",a2,"-",a2,"-",a2,
     .           ". Will treat this uby as a dummy. This term should not be sampled!")
        call printLine(WARNING,printstring)

c       now add this uby as a dummy to prevent future warnings
        iuby = NUByParam + 1
        
        if (iuby.le.MAXUBYPARAM) then
          tat1 = at1
          tat2 = at2 
          tat3 = at3
          call sortUreyBradley(tat1,tat2,tat3)
          
          UByName(iuby,1) = tat1
          UByName(iuby,2) = tat2
          UByName(iuby,3) = tat3
          UByParamIndex(iuby) = 0
          NUByParam = iuby
        endif

      endif

      if (findUreyBradleyParam.gt.0)
     .    findUreyBradleyParam = UByParamIndex(findUreyBradleyParam)
     
      return
      end
